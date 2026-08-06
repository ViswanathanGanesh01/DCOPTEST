import os
import shutil
import time
import warnings
import logging
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import fmpy
from fmpy import read_model_description, extract
from fmpy.util import write_csv
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing
from contextlib import contextmanager

# Suppress all Python and library warnings
warnings.filterwarnings('ignore')
logging.getLogger('fmpy').setLevel(logging.ERROR)

@contextmanager
def suppress_stdout_stderr():
    """A context manager that redirects stdout and stderr to devnull at the OS level."""
    devnull = os.open(os.devnull, os.O_WRONLY)
    saved_stdout_fd = os.dup(1)
    saved_stderr_fd = os.dup(2)
    try:
        os.dup2(devnull, 1)
        os.dup2(devnull, 2)
        yield
    finally:
        os.dup2(saved_stdout_fd, 1)
        os.dup2(saved_stderr_fd, 2)
        os.close(saved_stdout_fd)
        os.close(saved_stderr_fd)
        os.close(devnull)

# Global Settings
# Note: The cooling equipment in the FMUs is sized for a maximum IT load of 500 kW (500,000 W)
FMU_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "Model")
CLIMATES = ["1A"]
SIMULATION_STEP = 3600  # communication step size
DATA_STEP = 3600  #  DATA step size

# Day settings for simulation window
START_DAY = 0
END_DAY = 365

# Temperature Setpoints in Celsius (automatically converted to Kelvin in simulation)
CHW_SUP_SET = 15.0 # Chilled Water Supply Temp Setpoint (°C)
CDU_SUP_SET_C = 27.0   # CDU / Liquid Supply Temperature Setpoint (°C)
AIR_SUP_SET_C = 18.0   # Air Supply Temperature Setpoint (°C)
AIR_RET_SET_C = 25.0   # Air Return Temperature Setpoint (°C)
LIQ_RET_SET_C = 32.0   # Liquid Return Temperature Setpoint (°C)
CT_Range = 3    # Cooling Tower Range (°C)
T_FREEZE_C = 5.0       # Freeze Protection Temperature (°C)
T_APP_K = 6.0          # Approach Temperature Difference (K)

RESULTS_BASE_DIR = "results"
ERROR_LOG_PATH = "simulation_errors.log"
IT_LOAD_KW = 500.0  # IT load (kW) - Sized for max 500 kW

ARCHS_TO_RUN = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"]  # Explicitly configure architectures to simulate, evaluate, and plot
DRY_ARCH_IDS = ["4", "5", "6", "10", "11", "12"]
AIR_ARCH_IDS = [str(i) for i in range(1, 7)]
LIQUID_ARCH_IDS = [str(i) for i in range(7, 13)]

CARBON_INTENSITY = {"0A": 0.45, "Default": 0.40}

ARCH_NAME_MAP = {
    1: "Type 1", 2: "Type 2", 3: "Type 3", 4: "Type 4", 5: "Type 5", 6: "Type 6",
    7: "Type 7", 8: "Type 8", 9: "Type 9", 10: "Type 10", 11: "Type 11", 12: "Type 12"
}

ASHRAE_AIR_BINS = [(18, 27, "Rec"), (15, 32, "A1"), (10, 35, "A2"), (5, 40, "A3"), (5, 45, "A4")]
ASHRAE_LIQUID_BINS = [(2, 17, "W17"), (2, 27, "W27"), (2, 32, "W32"), (2, 40, "W40"), (2, 45, "W45")]
HATCHES = ['//', '\\\\', '||', '--', '++', 'xx', 'oo', '..', '**', 'OO']

# ----------------- Helper Functions -----------------

def get_enthalpy(temp_c, is_saturated=False, wb_c=None):
    def calc_psat(t):
        return 0.61121 * np.exp((18.678 - t / 234.5) * (t / (257.14 + t)))
    P_ATM = 101.325
    if is_saturated:
        p_sat = calc_psat(temp_c); w = 0.622 * p_sat / (P_ATM - p_sat)
        return 1.006 * temp_c + w * (2501 + 1.86 * temp_c), w
    if wb_c is not None:
        psat_wb = calc_psat(wb_c); w_s_wb = 0.622 * psat_wb / (P_ATM - psat_wb)
        h_fg_wb = 2501 - 2.326 * wb_c
        w = (h_fg_wb * w_s_wb - 1.006 * (temp_c - wb_c)) / (h_fg_wb + 1.86 * temp_c - 4.186 * wb_c)
        return 1.006 * temp_c + w * (2501 + 1.86 * temp_c), w
    return 0, 0

def calculate_pointwise_wue(df, it_load_kw):
    req = ["TCWRet", "TCWSup", "PumpCW", "Wetbulb", "Drybulb"]
    if not all(c in df.columns for c in req): return pd.Series([0.0] * len(df))
    COC, DRIFT_RATE, cp_water = 5.0, 0.0002, 4.186
    range_rr = df["TCWRet"] - df["TCWSup"]
    heat_rejected_kw = df["PumpCW"] * cp_water * range_rr
    wb_c = (df["Wetbulb"] - 273.15) if df["Wetbulb"].mean() > 200 else df["Wetbulb"]
    db_c = (df["Drybulb"] - 273.15) if df["Drybulb"].mean() > 200 else df["Drybulb"]
    h_in, w_in = get_enthalpy(db_c.values, wb_c=wb_c.values)
    t_exit = wb_c + 3.0
    h_out, w_out = get_enthalpy(t_exit.values, is_saturated=True)
    delta_h = np.clip(h_out - h_in, a_min=2.0, a_max=None)
    delta_w = np.clip(w_out - w_in, a_min=1e-5, a_max=None)
    evap = heat_rejected_kw * (delta_w / delta_h)
    makeup_kg_per_sec = evap + (evap / (COC - 1)) + (df["PumpCW"] * DRIFT_RATE)
    return (makeup_kg_per_sec * 3600) / it_load_kw

def get_compliance(series, bins):
    results = {}
    total = len(series)
    if total == 0: return {name: 0.0 for _, _, name in bins}
    for low, high, name in bins:
        count = ((series >= low) & (series <= high)).sum()
        results[name] = (count / total) * 100
    return results



def preprocess_climate_file(climate_file_path, processed_path, simulation_step):
    time_steps = np.arange(START_DAY * 24 * 3600.0, END_DAY * 24 * 3600.0 + simulation_step, simulation_step)
    df_climate = pd.read_csv(climate_file_path, sep='\t', header=None, comment='#', skiprows=lambda x: x == 1)
    climate_time = df_climate[0].values
    t_dry_c = df_climate[1].values
    rh = df_climate[3].values
    t_wet_c = t_dry_c * np.arctan(0.151977 * (rh + 8.313659)**0.5) + np.arctan(t_dry_c + rh) - np.arctan(rh - 1.676331) + 0.00391838 * (rh**1.5) * np.arctan(0.023101 * rh) - 4.686035
    
    tdry_profile = np.interp(time_steps, climate_time, t_dry_c) + 273.15
    twet_profile = np.interp(time_steps, climate_time, t_wet_c) + 273.15
    
    df_out = pd.DataFrame({
        'tdry_profile': tdry_profile,
        'twet_profile': twet_profile
    })
    df_out.to_csv(processed_path, index=False)

# ----------------- Parallel Worker -----------------

def run_single_simulation(args):
    arch_value, climate, current_sim_count, total_sims, processed_climate_path = args
    arch_sub_dir = f"Architecture_{int(arch_value)}"
    target_dir = os.path.join(RESULTS_BASE_DIR, arch_sub_dir)
    os.makedirs(target_dir, exist_ok=True)
    
    fmu_path = os.path.join(FMU_DIR, "DataCenterFMU.fmu")
    if not os.path.exists(fmu_path):
        return f"[{current_sim_count}/{total_sims}] Skipping: {climate} (DataCenterFMU.fmu not found)"
        
    climate_file_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "Climate Files", f"{climate}.mos")
    if not os.path.exists(climate_file_path):
        return f"[{current_sim_count}/{total_sims}] Skipping: {climate} (Climate file not found: {climate_file_path})"
        
    output_csv = os.path.join(target_dir, f"{climate}.csv")
    if os.path.exists(output_csv):
        return f"[{current_sim_count}/{total_sims}] Skipping: Arch {arch_value} | Climate {climate} (CSV exists)"

    extracted_fmu_dir = extract(fmu_path)
    try:
        model_description = read_model_description(extracted_fmu_dir, validate=False)
        prefixes = ["dataCenter.Temperature", "dataCenter.Control", "Temperature", "Control"]
        
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
            "control.Control.ValveChi": ["control.Control.ValveChi", "chiWSE.modulation.Chi.y"],
            "control.Control.ValveWSE": ["control.Control.ValveWSE", "chiWSE.modulation.WSE.y"]
        }
        
        model_var_names = [v.name for v in model_description.modelVariables]
        specific_outputs = []
        for std_name, name_options in candidates.items():
            found_var = None
            for opt in name_options:
                if opt in model_var_names:
                    found_var = opt
                    break
            if found_var: specific_outputs.append(found_var)
        
        output_vars = []
        for prefix in prefixes:
            output_vars.extend([v.name for v in model_description.modelVariables 
                                if v.name.startswith(f"{prefix}.") or v.name == prefix])
        output_vars.extend(specific_outputs)
        
        arch_var_name = 'architecture'
        if 'Arch' in model_var_names: arch_var_name = 'Arch'
        
        recorded_vars = ['time', arch_var_name]
        for v_name in ['TDryBulb', 'TWetBulb', 'ITAir', 'ITLiq']:
            if v_name in model_var_names: recorded_vars.append(v_name)
        recorded_vars.extend(output_vars)
        recorded_vars = list(dict.fromkeys(recorded_vars))

        model_inputs = [v.name for v in model_description.modelVariables if v.causality == 'input']
        all_possible_inputs = ['Arch', 'architecture', 'ITAir', 'ITLiq', 'AirSupSet', 'AirRetSet', 'LiqSupSet', 'LiquidRetSet', 'CWSupSet', 'CWRetSet', 'TFreeze', 'TApp', 'TCHWSupSet', 'TDryBulb', 'TWetBulb']
        input_dtype = [('time', np.float64)]
        for inp in all_possible_inputs:
            if inp in model_inputs: input_dtype.append((inp, np.float64))
        
        time_steps = np.arange(START_DAY * 24 * 3600.0, END_DAY * 24 * 3600.0 + SIMULATION_STEP, SIMULATION_STEP)
        inputs = np.zeros(len(time_steps), dtype=input_dtype)
        inputs['time'] = time_steps
        
        df_processed = pd.read_csv(processed_climate_path)
        tdry_profile = df_processed['tdry_profile'].values
        twet_profile = df_processed['twet_profile'].values
        
        if 'TDryBulb' in model_inputs: inputs['TDryBulb'] = tdry_profile
        if 'TWetBulb' in model_inputs: inputs['TWetBulb'] = twet_profile
        
        is_dry = str(int(arch_value)) in DRY_ARCH_IDS
        is_liquid = str(int(arch_value)) in LIQUID_ARCH_IDS
        chw_sup_set_c = CHW_SUP_SET
        
        if 'Arch' in model_inputs: inputs['Arch'] = arch_value
        elif 'architecture' in model_inputs: inputs['architecture'] = arch_value
            
        if 'AirSupSet' in model_inputs: inputs['AirSupSet'] = AIR_SUP_SET_C + 273.15
        if 'AirRetSet' in model_inputs: inputs['AirRetSet'] = AIR_RET_SET_C + 273.15
        if 'LiqSupSet' in model_inputs: inputs['LiqSupSet'] = CDU_SUP_SET_C + 273.15
        if 'LiquidRetSet' in model_inputs: inputs['LiquidRetSet'] = LIQ_RET_SET_C + 273.15
        if 'CWSupSet' in model_inputs:
            if is_dry:
                inputs['CWSupSet'] = np.minimum(23.89 + 273.15, np.maximum(12.78 + 273.15, tdry_profile + 3.0))
            else:
                inputs['CWSupSet'] = np.minimum(23.89 + 273.15, np.maximum(12.78 + 273.15, twet_profile + 3.0))
        if 'CWRetSet' in model_inputs:
            inputs['CWRetSet'] = inputs['CWSupSet'] + CT_Range
        if 'TFreeze' in model_inputs: inputs['TFreeze'] = T_FREEZE_C + 273.15
        if 'TApp' in model_inputs: inputs['TApp'] = 6.0 if is_dry else T_APP_K
        if 'TCHWSupSet' in model_inputs: inputs['TCHWSupSet'] = chw_sup_set_c + 273.15
            
        it_profile_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "it_profile.csv")
        custom_it_time = None
        if os.path.exists(it_profile_path):
            try:
                it_df = pd.read_csv(it_profile_path)
                custom_it_time = it_df['time'].values
                custom_it_air = it_df['ITAir'].values
                custom_it_liq = it_df['ITLiq'].values
            except Exception: pass
                
        if custom_it_time is not None:
            it_air_profile = np.interp(time_steps, custom_it_time, custom_it_air)
            it_liq_profile = np.interp(time_steps, custom_it_time, custom_it_liq)
        else:
            it_air_profile = np.full(len(time_steps), 120000.0)
            it_liq_profile = np.full(len(time_steps), 180000.0)
        
        if arch_value <= 6:
            it_air_final = it_air_profile
            it_liq_final = np.zeros(len(time_steps))
        else:
            it_air_final = np.zeros(len(time_steps))
            it_liq_final = it_liq_profile
        
        if 'ITAir' in model_inputs: inputs['ITAir'] = it_air_final
        if 'ITLiq' in model_inputs: inputs['ITLiq'] = it_liq_final

        start_vals = {}
        for var_name in ['Arch', 'architecture', 'AirSupSet', 'AirRetSet', 'LiqSupSet', 'LiquidRetSet', 'CWSupSet', 'CWRetSet', 'TFreeze', 'TApp', 'TCHWSupSet', 'TDryBulb', 'TWetBulb', 'ITAir', 'ITLiq', 'CDU.heatExchanger_TSet.k', 'CDU.heatExchanger_TSet.Ti', 'CDU.heatExchanger_TSet.Td']:
            if var_name in model_var_names:
                if var_name in ['Arch', 'architecture']: start_vals[var_name] = arch_value
                elif var_name == 'AirSupSet': start_vals[var_name] = AIR_SUP_SET_C + 273.15
                elif var_name == 'AirRetSet': start_vals[var_name] = AIR_RET_SET_C + 273.15
                elif var_name == 'LiqSupSet': start_vals[var_name] = CDU_SUP_SET_C + 273.15
                elif var_name == 'LiquidRetSet': start_vals[var_name] = LIQ_RET_SET_C + 273.15
                elif var_name == 'CWSupSet':
                    profile_val = tdry_profile[0] if is_dry else twet_profile[0]
                    min_lim = 12.78
                    start_vals[var_name] = min(23.89 + 273.15, max(min_lim + 273.15, profile_val + 3.0))
                elif var_name == 'CWRetSet':
                    start_vals[var_name] = start_vals['CWSupSet'] + 3
                elif var_name == 'TFreeze': start_vals[var_name] = T_FREEZE_C + 273.15
                elif var_name == 'TApp': start_vals[var_name] = 6.0 if is_dry else T_APP_K
                elif var_name == 'TCHWSupSet': start_vals[var_name] = chw_sup_set_c + 273.15
                elif var_name == 'TDryBulb': start_vals[var_name] = tdry_profile[0]
                elif var_name == 'TWetBulb': start_vals[var_name] = twet_profile[0]
                elif var_name == 'ITAir': start_vals[var_name] = it_air_final[0]
                elif var_name == 'ITLiq': start_vals[var_name] = it_liq_final[0]
                elif var_name == 'CDU.heatExchanger_TSet.k': start_vals[var_name] = 1.0
                elif var_name == 'CDU.heatExchanger_TSet.Ti': start_vals[var_name] = 60.0
                elif var_name == 'CDU.heatExchanger_TSet.Td': start_vals[var_name] = 10.0

        # Set tolerance to 0.1 to match Dymola's solver setting, resolving chattering/flow reversal instability
        dynamic_tol = 1e-4

        # Custom FMU logger to capture internal solver warnings/errors
        fmu_log_messages = []
        def fmu_logger(fmu, instanceName, status, category, message):
            msg = message.decode('utf-8', errors='ignore') if isinstance(message, bytes) else message
            if status >= 1:  # Status >= 1 indicates warning or error
                fmu_log_messages.append(msg)

        sim_kwargs = {
            'filename': extracted_fmu_dir, 'start_time': START_DAY * 24 * 3600.0, 'stop_time': END_DAY * 24 * 3600.0,
            'step_size': SIMULATION_STEP, 'output_interval': DATA_STEP,
            'model_description': model_description, 'output': recorded_vars,
            'relative_tolerance': dynamic_tol, 'logger': fmu_logger, 'start_values': start_vals
        }
        if len(input_dtype) > 1: sim_kwargs['input'] = inputs
        with suppress_stdout_stderr(): results = fmpy.simulate_fmu(**sim_kwargs)
        
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
            "dataCenter.pumCW.m_flow": "PumpCW", "pumCW.m_flow": "PumpCW", "control.Control.ValveChi": "ValveChi",
            "chiWSE.modulation.Chi.y": "ValveChi", "control.Control.ValveWSE": "ValveWSE", "chiWSE.modulation.WSE.y": "ValveWSE",
            "ITAir": "ITAir", "ITLiq": "ITLiq"
        }
        
        new_names = []
        seen_names = set()
        for name in results.dtype.names:
            cand = rename_map.get(name)
            if not cand:
                for p in prefixes:
                    if name.startswith(f"{p}."): cand = name[len(p)+1:]; break
            if not cand: cand = name
            if cand in ['Arch', 'architecture']: cand = 'architecture'
            final_name, counter = cand, 1
            while final_name in seen_names: final_name = f"{cand}_{counter}"; counter += 1
            new_names.append(final_name); seen_names.add(final_name)
        
        results.dtype.names = new_names
        df_res = pd.DataFrame(results)
        if 'PHVAC' not in df_res.columns:
            power_cols = ["Chiller_Power", "CoolingTower_FanPower", "PumpCW_Power", "PumpCHW_Power", "Fan_Power", "D2C_Pump_Power", "DCTFan_Power"]
            available_power_cols = [c for c in power_cols if c in df_res.columns]
            df_res['PHVAC'] = df_res[available_power_cols].sum(axis=1)
        if 'PUE' not in df_res.columns:
            it_air_val = df_res['ITAir'] if 'ITAir' in df_res.columns else 0.0
            it_liq_val = df_res['ITLiq'] if 'ITLiq' in df_res.columns else 0.0
            it_load = it_air_val + it_liq_val
            df_res['PUE'] = np.where(it_load > 0, 1.0 + df_res['PHVAC'] / it_load, 1.0)
            
        df_res.to_csv(output_csv, index=False)
        return f"[{current_sim_count}/{total_sims}] Success: Arch {arch_value} | Climate {climate}"
    except Exception as e:
        error_details = ""
        if fmu_log_messages:
            # Categorize common failure reasons from FMU C logs
            combined_logs = " | ".join(fmu_log_messages[-10:])
            reason = "Unknown solver issue"
            if "temperature" in combined_logs.lower() or "temp" in combined_logs.lower():
                if "exceeded" in combined_logs.lower() or "above" in combined_logs.lower() or "maximum" in combined_logs.lower():
                    reason = "Temperature above limit"
                elif "minimum" in combined_logs.lower() or "below" in combined_logs.lower():
                    reason = "Temperature below limit"
            elif "mass fraction" in combined_logs.lower() or "medium" in combined_logs.lower():
                reason = "Mass fraction limit bounds warning"
            elif "time step" in combined_logs.lower() or "stiff" in combined_logs.lower():
                reason = "Numerical stiffness / Step size too small"
            error_details = f" [Reason: {reason} | Last Logs: {combined_logs}]"
        
        error_msg = f"ERROR in Arch {arch_value}, Climate {climate}: {str(e)}{error_details}"
        with open(ERROR_LOG_PATH, "a") as f: f.write(error_msg + "\n")
        return f"[{current_sim_count}/{total_sims}] Failed: Arch {arch_value} | Climate {climate} - {str(e)}{error_details}"
    finally:
        shutil.rmtree(extracted_fmu_dir, ignore_errors=True)

def run_fmu_simulation():
    total_sims = len(ARCHS_TO_RUN) * len(CLIMATES)
    print(f"Initializing batch simulation: {total_sims} runs total.")
    print(f"Simulation step size set to: {SIMULATION_STEP} seconds.")
    
    processed_climate_paths = {}
    temp_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "Climate Files")
    for climate in CLIMATES:
        climate_file_path = os.path.join(temp_dir, f"{climate}.mos")
        processed_path = os.path.join(temp_dir, f"{climate}_temp_processed.csv")
        if os.path.exists(climate_file_path):
            print(f"[STATUS] Pre-processing climate data for {climate}...")
            preprocess_climate_file(climate_file_path, processed_path, SIMULATION_STEP)
            processed_climate_paths[climate] = processed_path

    try:
        tasks = []
        count = 0
        for arch_val_str in ARCHS_TO_RUN:
            arch_value = float(arch_val_str)
            for climate in CLIMATES:
                count += 1
                processed_path = processed_climate_paths.get(climate, "")
                tasks.append((arch_value, climate, count, total_sims, processed_path))
                
        num_cores = multiprocessing.cpu_count()
        print(f"Running simulations in parallel using {num_cores} CPU cores...")
        
        futures_map = {}
        completed_runs = 0
        with ProcessPoolExecutor(max_workers=num_cores) as executor:
            for t in tasks:
                arch_value, climate, count, total, _ = t
                arch_sub_dir = f"Architecture_{int(arch_value)}"
                output_csv = os.path.join(RESULTS_BASE_DIR, arch_sub_dir, f"{climate}.csv")
                if os.path.exists(output_csv):
                    print(f"[STATUS] Skipping (CSV already exists): Arch {arch_value} | Climate {climate} ({count}/{total})")
                    completed_runs += 1; continue
                print(f"[STATUS] Queueing: Arch {arch_value} | Climate {climate} ({count}/{total})")
                future = executor.submit(run_single_simulation, t)
                futures_map[future] = t
            if futures_map:
                for f in as_completed(futures_map):
                    arch_value, climate, count, total, _ = futures_map[f]
                    completed_runs += 1
                    try:
                        res = f.result()
                        print(f"[STATUS] ({completed_runs}/{total_sims}) Completed: {res}")
                    except Exception as e:
                        print(f"[STATUS] ({completed_runs}/{total_sims}) Exception in Arch {arch_value} | Climate {climate}: {e}")
            else:
                print(f"[STATUS] All {total_sims} simulations already completed (CSVs exist on disk).")
    finally:
        for path in processed_climate_paths.values():
            if os.path.exists(path):
                try: os.remove(path)
                except Exception: pass

def worker_helper(item):
    return extract_data_worker(item[0], item[1])

def extract_data_worker(arch_id, folder_name):
    results = []
    for climate in CLIMATES:
        file_path = os.path.join(RESULTS_BASE_DIR, f"Architecture_{arch_id}", f"{climate}.csv")
        if not os.path.exists(file_path): file_path = os.path.join(RESULTS_BASE_DIR, folder_name, f"{climate}.csv")
        if not os.path.exists(file_path): file_path = file_path.replace(".csv", ".CSV")
        if os.path.exists(file_path):
            try:
                df = pd.read_csv(file_path)
                pue = df['PUE'].mean() if 'PUE' in df.columns else 1.2
                wue = 0.0 if arch_id in DRY_ARCH_IDS else calculate_pointwise_wue(df, IT_LOAD_KW).mean()
                cue = pue * CARBON_INTENSITY.get(climate, CARBON_INTENSITY['Default'])
                fc_pct, pmc_pct, fmc_pct = 0, 0, 0
                if "ValveChi" in df.columns and "ValveWSE" in df.columns:
                    total_steps = len(df)
                    fc_pct = ((df["ValveChi"] == 0) & (df["ValveWSE"] == 1)).sum() / total_steps * 100
                    pmc_pct = ((df["ValveChi"] == 1) & (df["ValveWSE"] == 1)).sum() / total_steps * 100
                    fmc_pct = ((df["ValveChi"] == 1) & (df["ValveWSE"] == 0)).sum() / total_steps * 100
                temp_col = 'TAirSup' if arch_id in AIR_ARCH_IDS else 'TCDUSup'
                t_vals = df[temp_col] - 273.15 if temp_col in df.columns and df[temp_col].mean() > 200 else df.get(temp_col, pd.Series([]))
                
                # Initialize all bins to 0 first to prevent KeyError in plotting if only some archs are evaluated
                compliance = {}
                for _, _, b_name in ASHRAE_AIR_BINS:
                    compliance[b_name] = 0.0
                for _, _, b_name in ASHRAE_LIQUID_BINS:
                    compliance[b_name] = 0.0
                
                bins = ASHRAE_AIR_BINS if arch_id in AIR_ARCH_IDS else ASHRAE_LIQUID_BINS
                comp_vals = get_compliance(t_vals, bins)
                compliance.update(comp_vals)
                
                power_cols = {"Chiller": "Chiller_Power", "PumpCW": "PumpCW_Power", "PumpCHW": "PumpCHW_Power", "Fan": "Fan_Power", "D2CPump": "D2C_Pump_Power"}
                energy_vals = {label: (df[col].sum() * SIMULATION_STEP) / 3600000.0 if col in df.columns else 0.0 for label, col in power_cols.items()}
                
                # Combine both cooling tower and dry cooler fan power under a single CoolingTowerFan category
                ct_fan_energy = 0.0
                if "CoolingTower_FanPower" in df.columns:
                    ct_fan_energy += (df["CoolingTower_FanPower"].sum() * SIMULATION_STEP) / 3600000.0
                if "DCTFan_Power" in df.columns:
                    ct_fan_energy += (df["DCTFan_Power"].sum() * SIMULATION_STEP) / 3600000.0
                energy_vals["CoolingTowerFan"] = ct_fan_energy

                results.append({'ArchID': int(arch_id), 'ArchName': ARCH_NAME_MAP.get(int(arch_id)), 'Climate': climate, 'PUE': pue, 'WUE': wue, 'CUE': cue, 'FC': fc_pct, 'PMC': pmc_pct, 'FMC': fmc_pct, **compliance, **energy_vals})
            except Exception: pass
    return results

def plot_comprehensive_kpis(df):
    plt.style.use('seaborn-v0_8-white')
    fig, axes = plt.subplots(3, 2, figsize=(10, 7), sharex=True)
    plt.subplots_adjust(hspace=0.08, wspace=0.18)
    
    df = df.sort_values('ArchID')
    all_arch_ids = [int(x) for x in ARCHS_TO_RUN]
    x_labels = [ARCH_NAME_MAP[i] for i in all_arch_ids]

    kpi_metrics = [('PUE', 'PUE'), ('WUE', 'WUE (L/kWh)'), ('CUE', 'CUE (kgCO2/kWh)')]
    for i, (col, lab) in enumerate(kpi_metrics):
        ax = axes[i, 0]
        ax.grid(False)
        jitter = np.random.uniform(-0.15, 0.15, len(df))
        ax.scatter(df['ArchID'] - 1 + jitter, df[col], c=df[col], cmap='viridis', edgecolors='k', s=60, alpha=0.7, zorder=3)
        means = df.groupby('ArchID')[col].mean().reindex(all_arch_ids)
        ax.plot(np.arange(len(all_arch_ids)), means.values, 'r--', alpha=0.6, zorder=2)
        ax.set_ylabel(lab, fontweight='bold', fontsize=10)

    ax_m = axes[0, 1]
    ax_m.grid(False)
    mode_data = df.groupby('ArchID')[['FC', 'PMC', 'FMC']].mean().reindex(all_arch_ids)
    mode_data.plot(kind='bar', stacked=True, ax=ax_m, color=['#2ecc71', '#f1c40f', '#e74c3c'], width=0.75, edgecolor='black', linewidth=0.5, alpha=0.7)
    for i, container in enumerate(ax_m.containers):
        for bar in container:
            bar.set_hatch(HATCHES[i])
    ax_m.set_ylabel('Operating Mode (%)', fontweight='bold', fontsize=10)
    ax_m.set_ylim(0, 120)
    ax_m.legend(ncol=5, fontsize=9, loc='upper right')

    ax_air = axes[1, 1]
    ax_air.grid(False)
    air_cols = [b[2] for b in ASHRAE_AIR_BINS]
    comp_air = df.groupby('ArchID')[air_cols].mean().reindex(all_arch_ids)
    comp_air.plot(kind='bar', ax=ax_air, width=0.8, cmap='RdYlGn_r', edgecolor='black', linewidth=0.5, alpha=0.7)
    for i, container in enumerate(ax_air.containers):
        for bar in container:
            bar.set_hatch(HATCHES[i % len(HATCHES)])
    ax_air.set_ylabel('Air Compliance (%)', fontweight='bold', fontsize=10)
    ax_air.set_ylim(0, 115)
    ax_air.legend(ncol=5, fontsize=8, loc='upper right')

    ax_liq = axes[2, 1]
    ax_liq.grid(False)
    liq_cols = [b[2] for b in ASHRAE_LIQUID_BINS]
    comp_liq = df.groupby('ArchID')[liq_cols].mean().reindex(all_arch_ids)
    comp_liq.plot(kind='bar', ax=ax_liq, width=0.8, cmap='PuBuGn', edgecolor='black', linewidth=0.5, alpha=0.7)
    for i, container in enumerate(ax_liq.containers):
        for bar in container:
            bar.set_hatch(HATCHES[i % len(HATCHES)])
    ax_liq.set_ylabel('Liquid Compliance (%)', fontweight='bold', fontsize=10)
    ax_liq.set_ylim(0, 115)
    ax_liq.legend(ncol=4, fontsize=8, loc='upper left')

    for j in range(2):
        axes[2, j].set_xticks(range(len(all_arch_ids)))
        axes[2, j].set_xticklabels(x_labels, rotation=45, ha='right', fontweight='bold')
        axes[2, j].set_xlabel("Cooling System Type", fontweight='bold', fontsize=10)

    plt.savefig("KPI_Comprehensive_Analysis.pdf", dpi=300, bbox_inches='tight')

def plot_kpi_circle_analysis(df):
    plt.style.use('seaborn-v0_8-white')
    metrics = ['PUE', 'WUE', 'CUE']
    plot_archs = [int(x) for x in ARCHS_TO_RUN]
    summary_df = df.groupby('ArchID')[metrics].mean().reindex(plot_archs)
    
    max_vals = summary_df.max()
    min_vals = summary_df.min()
    
    norm_df = (summary_df - min_vals) / (max_vals - min_vals)
    angles = np.linspace(0, 2 * np.pi, len(metrics), endpoint=False)
    
    fig = plt.figure(figsize=(6, 4))
    ax = fig.add_subplot(111, polar=True)
    ax.grid(False)
    
    try:
        from matplotlib import colormaps
        colors = colormaps['tab20']
    except ImportError:
        colors = plt.cm.get_cmap('tab20', len(plot_archs))
    
    for i, arch_id in enumerate(plot_archs):
        vals = norm_df.loc[arch_id].values
        color = colors(i)
        
        x_coords = vals * np.cos(angles)
        y_coords = vals * np.sin(angles)
        centroid_r = np.sqrt(np.mean(x_coords)**2 + np.mean(y_coords)**2)
        centroid_theta = np.arctan2(np.mean(y_coords), np.mean(x_coords))
        
        plot_angles = np.append(angles, angles[0])
        plot_vals = np.append(vals, vals[0])
        ax.plot(plot_angles, plot_vals, color=color, linewidth=1.8, alpha=0.2, zorder=i+1, linestyle=":")
        ax.scatter(centroid_theta, centroid_r, color=color, s=100, edgecolors='black', 
                   alpha=1.0, zorder=50, label=ARCH_NAME_MAP[arch_id])

    ax.set_theta_offset(np.pi / 2)
    ax.set_theta_direction(-1)
    ax.set_xticks(angles)
    
    labels = []
    for m in metrics:
        mi, ma = min_vals[m], max_vals[m]
        labels.append(f"{m}\n(min:{mi:.2f}, max:{ma:.2f})")
            
    ax.set_xticklabels(labels, fontsize=11, fontweight='bold', ha='center')
    
    ax.tick_params(axis='x', which='major', pad=15) 
    for label in ax.get_xticklabels():
        if "WUE" in label.get_text() or "CUE" in label.get_text():
            label.set_y(-0.18) 

    ax.set_ylim(0, 1.1)
    ax.set_yticks([0.25, 0.5, 0.75, 1.0])
    ax.set_yticklabels(['25%', '50%', '75%', '100%'], fontsize=8, color='gray', fontweight='bold')

    for r in [0.25, 0.5, 0.75, 1.0]:
        circle = plt.Circle((0, 0), r, transform=ax.transData._b, color='gray', fill=False, linestyle=':', alpha=0.2)
        ax.add_artist(circle)
    for angle in angles:
        ax.plot([angle, angle], [0, 1], color='gray', linestyle=':', alpha=0.2, zorder=1)
    
    ax.spines['polar'].set_visible(True)
    ax.spines['polar'].set_linewidth(1.5)
    
    plt.legend(loc='upper right', bbox_to_anchor=(1.4, 1.05), title="Cooling Systems", fontsize=9)
    plt.savefig("KPI_Efficiency_Circle_Triangles.pdf", dpi=300, bbox_inches='tight')

def plot_equipment_energy(df):
    fig, ax = plt.subplots(figsize=(10, 6))
    energy_cols = ["Chiller", "CoolingTowerFan", "PumpCW", "PumpCHW", "Fan", "D2CPump"]
    plot_archs = [int(x) for x in ARCHS_TO_RUN]
    data = df.groupby('ArchID')[energy_cols].mean().reindex(plot_archs) / 1000.0
    data.rename(columns={"Chiller": "Chiller", "CoolingTowerFan": "Cooling Tower Fan", "PumpCW": "Condenser Water Pump", "PumpCHW": "Chilled Water Pump", "Fan": "AHU Fan", "D2CPump": "D2C Pump"}).plot(kind='bar', stacked=True, ax=ax, width=0.75, edgecolor='black', alpha=0.8)
    for i, container in enumerate(ax.containers):
        for bar in container:
            bar.set_hatch(HATCHES[i % len(HATCHES)])
    ax.set_ylabel('Energy (MWh)', fontweight='bold')
    x_labels = [ARCH_NAME_MAP[i] for i in plot_archs]
    ax.set_xticks(range(len(plot_archs)))
    ax.set_xticklabels(x_labels, rotation=45, ha='right', fontweight='bold')
    ax.set_xlabel('Cooling System Type', fontweight='bold', fontsize=10)
    max_val = data.sum(axis=1).max()
    ax.set_ylim(0, 1.1 * max_val)
    plt.savefig("KPI_Equipment_Energy_Consumption.pdf", dpi=300, bbox_inches='tight')
    

def main():
    try:
        run_fmu_simulation()
        folders = [f for f in os.listdir(RESULTS_BASE_DIR) if os.path.isdir(os.path.join(RESULTS_BASE_DIR, f))]
        plot_archs = [int(x) for x in ARCHS_TO_RUN]
        all_map = { "".join(filter(str.isdigit, f)): f for f in folders if int("".join(filter(str.isdigit, f))) in plot_archs }
        all_data = []
        with ProcessPoolExecutor(max_workers=multiprocessing.cpu_count()) as executor:
            for res in executor.map(worker_helper, all_map.items()):
                if res: all_data.extend(res)
        if all_data:
            df = pd.DataFrame(all_data)
            plot_comprehensive_kpis(df); plot_kpi_circle_analysis(df); plot_equipment_energy(df)
    finally:
        pass

if __name__ == "__main__":
    multiprocessing.freeze_support()
    main()