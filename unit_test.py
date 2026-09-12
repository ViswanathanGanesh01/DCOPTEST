"""
Unit test suite for DataCenterFMU and cooling simulation pipeline.
Compatible with both Windows and Linux environments.

Includes:
- Environment & package dependency checks
- FMU archive and platform binary integrity checks (.dll on Windows, .so on Linux)
- Psychrometric and WUE calculation checks
- Simulation smoke execution (24 hours)
- Numerical validation against reference baseline results with configurable tolerances

Run with:
    python test_simulation.py
or
    python -m unittest test_simulation.py
"""

import os
import sys
import shutil
import zipfile
import unittest
import numpy as np
import pandas as pd
import fmpy
from fmpy import read_model_description, extract

import warnings
import logging

# Suppress all Python, FMPy, and C-level warnings
warnings.filterwarnings('ignore')
logging.getLogger('fmpy').setLevel(logging.CRITICAL)
logging.disable(logging.CRITICAL)

# Add BASE_DIR to sys.path
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
if BASE_DIR not in sys.path:
    sys.path.insert(0, BASE_DIR)

import Example

MODEL_DIR = os.path.join(BASE_DIR, "Model")
CLIMATE_DIR = os.path.join(BASE_DIR, "Climate Files")
TEST_DATA_DIR = os.path.join(BASE_DIR, "test_data")
FMU_PATH = os.path.join(MODEL_DIR, "DataCenterFMU.fmu")
IT_PROFILE_PATH = os.path.join(BASE_DIR, "it_profile.csv")
REF_DATA_PATH = os.path.join(TEST_DATA_DIR, "ref_arch1_24h.csv")


def run_smoke_simulation_df(arch_value=1.0, hours=24):
    """
    Executes a short simulation and returns the resulting DataFrame
    with standard column names (PHVAC, PUE, powers, etc.).
    """
    climate_mos = os.path.join(CLIMATE_DIR, "1A.mos")
    step_size = 3600.0
    time_steps = np.arange(0, hours * 3600.0 + step_size, step_size)

    df_climate = pd.read_csv(climate_mos, sep='\t', header=None, comment='#', skiprows=lambda x: x == 1)
    c_time = df_climate[0].values
    t_dry = df_climate[1].values
    rh = df_climate[3].values
    t_wet = t_dry * np.arctan(0.151977 * (rh + 8.313659)**0.5) + np.arctan(t_dry + rh) - np.arctan(rh - 1.676331) + 0.00391838 * (rh**1.5) * np.arctan(0.023101 * rh) - 4.686035

    tdry_profile = np.interp(time_steps, c_time, t_dry) + 273.15
    twet_profile = np.interp(time_steps, c_time, t_wet) + 273.15

    extracted_dir = extract(FMU_PATH)
    try:
        model_description = read_model_description(extracted_dir, validate=False)
        model_var_names = [v.name for v in model_description.modelVariables]
        model_inputs = [v.name for v in model_description.modelVariables if v.causality == 'input']

        all_possible_inputs = ['Arch', 'architecture', 'ITAir', 'ITLiq', 'AirSupSet', 'AirRetSet', 'LiqSupSet', 'LiquidRetSet', 'CWSupSet', 'CWRetSet', 'TFreeze', 'TApp', 'TCHWSupSet', 'TDryBulb', 'TWetBulb']
        input_dtype = [('time', np.float64)]
        for inp in all_possible_inputs:
            if inp in model_inputs:
                input_dtype.append((inp, np.float64))

        inputs = np.zeros(len(time_steps), dtype=input_dtype)
        inputs['time'] = time_steps

        if 'TDryBulb' in model_inputs: inputs['TDryBulb'] = tdry_profile
        if 'TWetBulb' in model_inputs: inputs['TWetBulb'] = twet_profile
        if 'Arch' in model_inputs: inputs['Arch'] = arch_value
        if 'architecture' in model_inputs: inputs['architecture'] = arch_value
        if 'ITAir' in model_inputs: inputs['ITAir'] = 120000.0
        if 'ITLiq' in model_inputs: inputs['ITLiq'] = 0.0
        if 'AirSupSet' in model_inputs: inputs['AirSupSet'] = 18.0 + 273.15
        if 'AirRetSet' in model_inputs: inputs['AirRetSet'] = 25.0 + 273.15
        if 'LiqSupSet' in model_inputs: inputs['LiqSupSet'] = 27.0 + 273.15
        if 'LiquidRetSet' in model_inputs: inputs['LiquidRetSet'] = 32.0 + 273.15
        if 'CWSupSet' in model_inputs:
            inputs['CWSupSet'] = np.minimum(23.89 + 273.15, np.maximum(12.78 + 273.15, twet_profile + 3.0))
        if 'CWRetSet' in model_inputs:
            inputs['CWRetSet'] = inputs['CWSupSet'] + 3.0
        if 'TFreeze' in model_inputs: inputs['TFreeze'] = 5.0 + 273.15
        if 'TApp' in model_inputs: inputs['TApp'] = 6.0
        if 'TCHWSupSet' in model_inputs: inputs['TCHWSupSet'] = 15.0 + 273.15

        start_vals = {}
        for var_name in ['Arch', 'architecture', 'AirSupSet', 'AirRetSet', 'LiqSupSet', 'LiquidRetSet', 'CWSupSet', 'CWRetSet', 'TFreeze', 'TApp', 'TCHWSupSet', 'TDryBulb', 'TWetBulb', 'ITAir', 'ITLiq', 'CDU.heatExchanger_TSet.k', 'CDU.heatExchanger_TSet.Ti', 'CDU.heatExchanger_TSet.Td']:
            if var_name in model_var_names:
                if var_name in ['Arch', 'architecture']: start_vals[var_name] = arch_value
                elif var_name == 'AirSupSet': start_vals[var_name] = 18.0 + 273.15
                elif var_name == 'AirRetSet': start_vals[var_name] = 25.0 + 273.15
                elif var_name == 'LiqSupSet': start_vals[var_name] = 27.0 + 273.15
                elif var_name == 'LiquidRetSet': start_vals[var_name] = 32.0 + 273.15
                elif var_name == 'CWSupSet': start_vals[var_name] = float(inputs['CWSupSet'][0]) if 'CWSupSet' in model_inputs else 20.0 + 273.15
                elif var_name == 'CWRetSet': start_vals[var_name] = start_vals.get('CWSupSet', 20.0 + 273.15) + 3.0
                elif var_name == 'TFreeze': start_vals[var_name] = 5.0 + 273.15
                elif var_name == 'TApp': start_vals[var_name] = 6.0
                elif var_name == 'TCHWSupSet': start_vals[var_name] = 15.0 + 273.15
                elif var_name == 'TDryBulb': start_vals[var_name] = float(tdry_profile[0])
                elif var_name == 'TWetBulb': start_vals[var_name] = float(twet_profile[0])
                elif var_name == 'ITAir': start_vals[var_name] = 120000.0
                elif var_name == 'ITLiq': start_vals[var_name] = 0.0

        candidates = {
            "dataCenter.PHVAC.y": ["dataCenter.PHVAC.y"],
            "dataCenter.PUE.y": ["dataCenter.PUE.y"],
            "dataCenter.pumCW.P": ["dataCenter.pumCW.P", "pumCW.P"],
            "dataCenter.DCTFan.P": ["dataCenter.DCTFan.P", "DCTFan.P"],
            "dataCenter.cooTow.PFan": ["dataCenter.cooTow.PFan", "cooTow.PFan"],
            "dataCenter.pumCHW.P": ["dataCenter.pumCHW.P", "integratedPrimaryLoadSide.pum.pum[1].P"],
            "dataCenter.fan.P": ["dataCenter.ahu.fan.P", "dataCenter.fan.P", "AHU.fan.P"],
            "dataCenter.d2C.pum.P": ["dataCenter.DirectToChip.pum.P", "dataCenter.d2C.pum.P", "CDU.pum.P"],
            "dataCenter.chi.P": ["dataCenter.chi.P", "integratedPrimaryLoadSide.chiPar.chi[1].P"],
            "dataCenter.d2C.TRet.T": ["dataCenter.DirectToChip.TRet.T", "dataCenter.d2C.TRet.T", "Temperature.TCDURet"],
            "control.weaBus.TDryBul": ["control.weaBus.TDryBul", "TDryBulb"],
            "control.weaBus.TWetBul": ["control.weaBus.TWetBul", "TWetBulb"],
            "dataCenter.pumCW.m_flow": ["dataCenter.pumCW.m_flow", "pumCW.m_flow"],
        }
        
        specific_outputs = []
        for std_name, name_options in candidates.items():
            for opt in name_options:
                if opt in model_var_names:
                    specific_outputs.append(opt)
                    break

        recorded_vars = list(dict.fromkeys(['time', 'architecture'] + specific_outputs))

        def silent_fmu_logger(fmu, instanceName, status, category, message):
            pass

        with Example.suppress_stdout_stderr():
            results = fmpy.simulate_fmu(
                filename=extracted_dir,
                start_time=0.0,
                stop_time=hours * 3600.0,
                step_size=step_size,
                output_interval=step_size,
                input=inputs,
                output=recorded_vars,
                start_values=start_vals,
                relative_tolerance=1e-4,
                logger=silent_fmu_logger
            )

        rename_map = {
            "dataCenter.PHVAC.y": "PHVAC", "dataCenter.PUE.y": "PUE",
            "dataCenter.pumCW.P": "PumpCW_Power", "pumCW.P": "PumpCW_Power",
            "dataCenter.DCTFan.P": "DCTFan_Power", "DCTFan.P": "DCTFan_Power",
            "dataCenter.cooTow.PFan": "CoolingTower_FanPower", "cooTow.PFan": "CoolingTower_FanPower",
            "dataCenter.pumCHW.P": "PumpCHW_Power", "integratedPrimaryLoadSide.pum.pum[1].P": "PumpCHW_Power",
            "dataCenter.fan.P": "Fan_Power", "dataCenter.ahu.fan.P": "Fan_Power", "AHU.fan.P": "Fan_Power",
            "dataCenter.d2C.pum.P": "D2C_Pump_Power", "dataCenter.DirectToChip.pum.P": "D2C_Pump_Power", "CDU.pum.P": "D2C_Pump_Power",
            "dataCenter.chi.P": "Chiller_Power", "integratedPrimaryLoadSide.chiPar.chi[1].P": "Chiller_Power",
            "dataCenter.d2C.TRet.T": "TCDURet", "dataCenter.DirectToChip.TRet.T": "TCDURet", "Temperature.TCDURet": "TCDURet",
            "control.weaBus.TDryBul": "Drybulb", "TDryBulb": "Drybulb", "control.weaBus.TWetBul": "Wetbulb", "TWetBulb": "Wetbulb",
            "dataCenter.pumCW.m_flow": "PumpCW", "pumCW.m_flow": "PumpCW",
            "ITAir": "ITAir", "ITLiq": "ITLiq"
        }

        new_names = [rename_map.get(n, n) for n in results.dtype.names]
        results.dtype.names = new_names
        df_res = pd.DataFrame(results)

        if 'PHVAC' not in df_res.columns:
            power_cols = ["Chiller_Power", "CoolingTower_FanPower", "PumpCW_Power", "PumpCHW_Power", "Fan_Power", "D2C_Pump_Power", "DCTFan_Power"]
            available_power_cols = [c for c in power_cols if c in df_res.columns]
            df_res['PHVAC'] = df_res[available_power_cols].sum(axis=1)

        if 'PUE' not in df_res.columns:
            it_load = 120000.0
            df_res['PUE'] = 1.0 + df_res['PHVAC'] / it_load

        return df_res

    finally:
        shutil.rmtree(extracted_dir, ignore_errors=True)


class TestEnvironmentAndDependencies(unittest.TestCase):
    """Verifies environment setup and required packages."""

    def test_dependencies(self):
        for pkg in ["numpy", "pandas", "matplotlib", "fmpy"]:
            with self.subTest(package=pkg):
                __import__(pkg)

    def test_workspace_files_exist(self):
        self.assertTrue(os.path.isdir(MODEL_DIR), "Model directory missing")
        self.assertTrue(os.path.isdir(CLIMATE_DIR), "Climate Files directory missing")
        self.assertTrue(os.path.isfile(FMU_PATH), f"FMU file missing: {FMU_PATH}")
        self.assertTrue(os.path.isfile(IT_PROFILE_PATH), f"IT profile missing: {IT_PROFILE_PATH}")


class TestFMUIntegrity(unittest.TestCase):
    """Validates FMU archive structure, XML metadata, and platform binaries."""

    def setUp(self):
        self.assertTrue(os.path.exists(FMU_PATH), "FMU file not found")
        self.zip_file = zipfile.ZipFile(FMU_PATH)

    def tearDown(self):
        self.zip_file.close()

    def test_model_description_xml_present(self):
        namelist = self.zip_file.namelist()
        self.assertIn("modelDescription.xml", namelist, "modelDescription.xml missing from FMU")

    def test_platform_binary_present(self):
        namelist = self.zip_file.namelist()
        if sys.platform == "win32":
            has_binary = ("binaries/win64/DataCenterFMU.dll" in namelist or 
                          "binaries/win32/DataCenterFMU.dll" in namelist)
            self.assertTrue(has_binary, "Windows DLL missing in FMU binaries/")
        elif sys.platform.startswith("linux"):
            has_binary = any(n.startswith("binaries/linux64/") and n.endswith(".so") for n in namelist)
            self.assertTrue(has_binary, "Linux .so binary missing in FMU binaries/linux64/")

    def test_read_model_description(self):
        temp_extracted = extract(FMU_PATH)
        try:
            desc = read_model_description(temp_extracted, validate=False)
            self.assertEqual(desc.fmiVersion, "2.0", "FMU must be FMI 2.0")
            var_names = [v.name for v in desc.modelVariables]
            
            # Check for essential input/architecture variables
            self.assertTrue("Arch" in var_names or "architecture" in var_names, 
                            "Architecture input variable missing")
            for expected in ["TDryBulb", "TWetBulb", "ITAir", "ITLiq"]:
                self.assertIn(expected, var_names, f"Expected variable '{expected}' missing")
        finally:
            shutil.rmtree(temp_extracted, ignore_errors=True)


class TestPsychrometricCalculations(unittest.TestCase):
    """Tests enthalpy and psychrometric functions from Example.py."""

    def setUp(self):
        from Example import get_enthalpy, calculate_pointwise_wue
        self.get_enthalpy = get_enthalpy
        self.calculate_pointwise_wue = calculate_pointwise_wue

    def test_enthalpy_reasonable_range(self):
        # At 25°C dry bulb and 18°C wet bulb, enthalpy should be approx 50-60 kJ/kg
        h, w = self.get_enthalpy(25.0, wb_c=18.0)
        self.assertGreater(h, 40.0)
        self.assertLess(h, 70.0)
        self.assertGreater(w, 0.005)
        self.assertLess(w, 0.020)

    def test_wue_calculation(self):
        df_dummy = pd.DataFrame({
            "TCWRet": [30.0, 32.0],
            "TCWSup": [25.0, 26.0],
            "PumpCW": [10.0, 10.0],
            "Wetbulb": [20.0 + 273.15, 22.0 + 273.15],
            "Drybulb": [30.0 + 273.15, 32.0 + 273.15]
        })
        wue_series = self.calculate_pointwise_wue(df_dummy, it_load_kw=500.0)
        self.assertEqual(len(wue_series), 2)
        self.assertTrue((wue_series >= 0).all(), "WUE must be non-negative")


class TestSimulationSmoke(unittest.TestCase):
    """Performs a fast 24-hour smoke simulation of Architecture 1."""

    def test_24hr_simulation(self):
        df = run_smoke_simulation_df(arch_value=1.0, hours=24)
        self.assertIsNotNone(df, "Simulation DataFrame should not be None")
        self.assertEqual(len(df), 25, "Expected 25 hourly data points for 24h")
        self.assertIn("PHVAC", df.columns, "PHVAC column missing")
        self.assertIn("PUE", df.columns, "PUE column missing")
        
        # Verify no NaN or Inf in outputs
        self.assertFalse(df.isna().any().any(), "Output contains NaN values")
        self.assertFalse(np.isinf(df.to_numpy()).any(), "Output contains Inf values")


class TestReferenceComparison(unittest.TestCase):
    """Compares fresh simulation output against stored golden reference results."""

    @classmethod
    def setUpClass(cls):
        if not os.path.isfile(REF_DATA_PATH):
            print(f"\n[REFERENCE CHECK] Reference data file not found at: {REF_DATA_PATH}")
            should_write = False
            if sys.stdin:
                try:
                    ans = input("No reference data exists. Would you like to run the simulation and write it as test_data? [y/N]: ").strip().lower()
                    should_write = ans in ("y", "yes")
                except (EOFError, KeyboardInterrupt):
                    should_write = False
            elif os.environ.get("AUTO_WRITE_REF_DATA", "").lower() in ("1", "true", "yes"):
                should_write = True

            if should_write:
                print("[REFERENCE CHECK] Running simulation to generate reference test_data...")
                os.makedirs(os.path.dirname(REF_DATA_PATH), exist_ok=True)
                new_ref_df = run_smoke_simulation_df(arch_value=1.0, hours=24)
                new_ref_df.to_csv(REF_DATA_PATH, index=False)
                print(f"[REFERENCE CHECK] Golden reference successfully saved to: {REF_DATA_PATH}")
            else:
                raise unittest.SkipTest(
                    f"Reference dataset missing at '{REF_DATA_PATH}' and generation was skipped or declined."
                )

        cls.ref_df = pd.read_csv(REF_DATA_PATH)

    def test_compare_simulation_against_reference(self):
        sim_df = run_smoke_simulation_df(arch_value=1.0, hours=24)
        
        # 1. Verify shape match
        self.assertEqual(len(sim_df), len(self.ref_df), "Row count differs from reference dataset")
        
        # 2. Tolerances for physical quantities
        tolerances = {
            "PHVAC":                 {"rtol": 1e-3, "atol": 5.0},    # Cooling power (W)
            "PUE":                   {"rtol": 1e-4, "atol": 0.005},  # Dimensionless PUE
            "Fan_Power":             {"rtol": 1e-3, "atol": 2.0},    # AHU Fan Power (W)
            "Chiller_Power":         {"rtol": 1e-3, "atol": 2.0},    # Chiller Power (W)
            "PumpCW_Power":          {"rtol": 1e-3, "atol": 1.0},    # Condenser Pump Power (W)
            "CoolingTower_FanPower": {"rtol": 1e-3, "atol": 1.0},    # CT Fan Power (W)
            "Drybulb":               {"rtol": 1e-4, "atol": 0.05},   # Ambient Drybulb (K)
            "Wetbulb":               {"rtol": 1e-4, "atol": 0.05},   # Ambient Wetbulb (K)
        }

        for metric, tol in tolerances.items():
            if metric in sim_df.columns and metric in self.ref_df.columns:
                with self.subTest(metric=metric):
                    np.testing.assert_allclose(
                        sim_df[metric].to_numpy(),
                        self.ref_df[metric].to_numpy(),
                        rtol=tol["rtol"],
                        atol=tol["atol"],
                        err_msg=f"Metric '{metric}' diverged beyond allowable tolerance from reference!"
                    )


if __name__ == "__main__":
    unittest.main(verbosity=2)
