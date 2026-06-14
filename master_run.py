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
# Note: The cooling equipment in the FMUs is sized for a maximum IT load of 1 MW (1,000,000 W)
# Note: Currently runs for a single climate. Additional climates will be explored in the future.
FMU_DIR = os.path.dirname(os.path.abspath(__file__))
CLIMATES = ["DataCenter_0A"]
SIMULATION_STEP = 600  # 10 minute communication step size

# Temperature Setpoints in Celsius (automatically converted to Kelvin in simulation)
CHW_SUP_SET_C = 15.0  # Chilled Water Supply Temperature Setpoint (°C)
CDU_SUP_SET_C = 27.0  # CDU Supply Temperature Setpoint (°C)

RESULTS_BASE_DIR = "results"
ERROR_LOG_PATH = "simulation_errors.log"
IT_LOAD_KW = 180.0  # IT load (kW) - Sized for max 1000 kW (1 MW)

ALL_ARCHS = [str(i) for i in range(1, 13)]
DRY_ARCH_IDS = ["4", "5", "6", "10", "11", "12"]
AIR_ARCH_IDS = [str(i) for i in range(1, 7)]
LIQUID_ARCH_IDS = [str(i) for i in range(7, 13)]

CARBON_INTENSITY = {"DataCenter_0A": 0.45, "DataCenter_0B": 0.45, "Default": 0.40}

ARCH_NAME_MAP = {
    1: "Type 1",
    2: "Type 2",
    3: "Type 3",
    4: "Type 4",
    5: "Type 5",
    6: "Type 6",
    7: "Type 7",
    8: "Type 8",
    9: "Type 9",
    10: "Type 10",
    11: "Type 11",
    12: "Type 12"
}

ASHRAE_AIR_BINS = [
    (18, 27, "Rec"), (15, 32, "A1"), (10, 35, "A2"), (5, 40, "A3"), (5, 45, "A4")
]
ASHRAE_LIQUID_BINS = [
    (2, 17, "W17"), (2, 27, "W27"), (2, 32, "W32"), (2, 40, "W40"), (2, 45, "W45")
]

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

# ----------------- Parallel Worker -----------------

def run_single_simulation(args):
    arch_value, climate, current_sim_count, total_sims = args
    arch_sub_dir = f"Architecture_{int(arch_value)}"
    target_dir = os.path.join(RESULTS_BASE_DIR, arch_sub_dir)
    os.makedirs(target_dir, exist_ok=True)
    
    fmu_path = os.path.join(FMU_DIR, f"{climate}.fmu")
    if not os.path.exists(fmu_path):
        return f"[{current_sim_count}/{total_sims}] Skipping: {climate} (File not found)"
        
    output_csv = os.path.join(target_dir, f"{climate}.csv")
    if os.path.exists(output_csv):
        return f"[{current_sim_count}/{total_sims}] Skipping: Arch {arch_value} | Climate {climate} (CSV exists)"

    unzipdir = None
    try:
        model_description = read_model_description(fmu_path)
        prefixes = ["dataCenter.Temperature", "dataCenter.Control"]
        
        candidates = {
            "dataCenter.PHVAC.y": ["dataCenter.PHVAC.y"],
            "dataCenter.PUE.y": ["dataCenter.PUE.y"],
            "dataCenter.pumCW.P": ["dataCenter.pumCW.P"],
            "dataCenter.DCTFan.P": ["dataCenter.DCTFan.P"],
            "dataCenter.cooTow.PFan": ["dataCenter.cooTow.PFan"],
            "dataCenter.pumCHW.P": ["dataCenter.pumCHW.P"],
            "dataCenter.fan.P": ["dataCenter.ahu.fan.P", "dataCenter.fan.P"],
            "dataCenter.d2C.pum.P": ["dataCenter.DirectToChip.pum.P", "dataCenter.d2C.pum.P"],
            "dataCenter.chi.P": ["dataCenter.chi.P"],
            "dataCenter.d2C.TRet.T": ["dataCenter.DirectToChip.TRet.T", "dataCenter.d2C.TRet.T"],
            "control.weaBus.TDryBul": ["control.weaBus.TDryBul"],
            "control.weaBus.TWetBul": ["control.weaBus.TWetBul"],
            "dataCenter.pumCW.m_flow": ["dataCenter.pumCW.m_flow"],
            "control.Control.ValveChi": ["control.Control.ValveChi"],
            "control.Control.ValveWSE": ["control.Control.ValveWSE"]
        }
        
        model_var_names = [v.name for v in model_description.modelVariables]
        
        specific_outputs = []
        for std_name, name_options in candidates.items():
            found_var = None
            for opt in name_options:
                if opt in model_var_names:
                    found_var = opt
                    break
            if found_var:
                specific_outputs.append(found_var)
        
        output_vars = []
        for prefix in prefixes:
            output_vars.extend([v.name for v in model_description.modelVariables 
                              if v.name.startswith(f"{prefix}.") or v.name == prefix])
        output_vars.extend(specific_outputs)
        
        arch_var_name = 'architecture'
        if 'Arch' in model_var_names:
            arch_var_name = 'Arch'
        
        recorded_vars = ['time', arch_var_name] + output_vars

        # Inputs
        model_inputs = [v.name for v in model_description.modelVariables if v.causality == 'input']
        
        input_dtype = [('time', np.float64)]
        for inp in ['Arch', 'architecture', 'CHWSupSet', 'CDUSupSet', 'ITAir', 'ITLiq']:
            if inp in model_inputs:
                input_dtype.append((inp, np.float64))
        
        time_steps = np.arange(0.0, 365 * 24 * 3600 + SIMULATION_STEP, SIMULATION_STEP)
        inputs = np.zeros(len(time_steps), dtype=input_dtype)
        inputs['time'] = time_steps
        
        if 'Arch' in model_inputs:
            inputs['Arch'] = arch_value
        elif 'architecture' in model_inputs:
            inputs['architecture'] = arch_value
            
        if 'CHWSupSet' in model_inputs:
            inputs['CHWSupSet'] = CHW_SUP_SET_C + 273.15
        if 'CDUSupSet' in model_inputs:
            inputs['CDUSupSet'] = CDU_SUP_SET_C + 273.15
            
        # Load custom IT profile
        it_profile_path = os.path.join(FMU_DIR, "it_profile.csv")
        custom_it_time = None
        if os.path.exists(it_profile_path):
            try:
                it_df = pd.read_csv(it_profile_path)
                custom_it_time = it_df['time'].values
                custom_it_air = it_df['ITAir'].values
                custom_it_liq = it_df['ITLiq'].values
            except Exception:
                custom_it_time = None
                
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
        
        if 'ITAir' in model_inputs:
            inputs['ITAir'] = it_air_final
        if 'ITLiq' in model_inputs:
            inputs['ITLiq'] = it_liq_final

        unzipdir = extract(fmu_path)
        
        sim_kwargs = {
            'filename': unzipdir,
            'start_time': 0.0,
            'stop_time': 365 * 24 * 3600,
            'step_size': SIMULATION_STEP,
            'output_interval': SIMULATION_STEP,
            'model_description': model_description,
            'output': recorded_vars,
            'relative_tolerance': 1e-5,
            'logger': lambda *args: None
        }
        
        if len(input_dtype) > 1:
            sim_kwargs['input'] = inputs
        else:
            sim_kwargs['start_values'] = {arch_var_name: arch_value}
        
        # Execute under OS-level output suppression to clean Newton warnings
        with suppress_stdout_stderr():
            results = fmpy.simulate_fmu(**sim_kwargs)
        
        rename_map = {
            "dataCenter.PHVAC.y": "PHVAC", 
            "dataCenter.PUE.y": "PUE",
            "dataCenter.pumCW.P": "PumpCW_Power",
            "dataCenter.DCTFan.P": "DCTFan_Power",
            "dataCenter.cooTow.PFan": "CoolingTower_FanPower",
            "dataCenter.pumCHW.P": "PumpCHW_Power",
            "dataCenter.fan.P": "Fan_Power",
            "dataCenter.ahu.fan.P": "Fan_Power",
            "dataCenter.d2C.pum.P": "D2C_Pump_Power",
            "dataCenter.DirectToChip.pum.P": "D2C_Pump_Power",
            "dataCenter.chi.P": "Chiller_Power",
            "dataCenter.d2C.TRet.T": "TCDURet",
            "dataCenter.DirectToChip.TRet.T": "TCDURet",
            "control.weaBus.TDryBul": "Drybulb",
            "control.weaBus.TWetBul": "Wetbulb",
            "dataCenter.pumCW.m_flow": "PumpCW",
            "control.Control.ValveChi": "ValveChi",
            "control.Control.ValveWSE": "ValveWSE"
        }
        
        new_names = []
        seen_names = set()
        
        for name in results.dtype.names:
            cand = rename_map.get(name)
            if not cand:
                for p in prefixes:
                    if name.startswith(f"{p}."):
                        cand = name[len(p)+1:]
                        break
            if not cand: cand = name
            if cand in ['Arch', 'architecture']:
                cand = 'architecture'
            
            final_name, counter = cand, 1
            while final_name in seen_names:
                final_name = f"{cand}_{counter}"
                counter += 1
            new_names.append(final_name)
            seen_names.add(final_name)
        
        results.dtype.names = new_names
        write_csv(output_csv, results)
        return f"[{current_sim_count}/{total_sims}] Success: Arch {arch_value} | Climate {climate}"
        
    except Exception as e:
        error_msg = f"ERROR in Arch {arch_value}, Climate {climate}: {str(e)}"
        with open(ERROR_LOG_PATH, "a") as f:
            f.write(error_msg + "\n")
        return f"[{current_sim_count}/{total_sims}] Failed: Arch {arch_value} | Climate {climate} - {str(e)}"
    
    finally:
        if unzipdir and os.path.exists(unzipdir):
            shutil.rmtree(unzipdir)

def run_fmu_simulation():
    total_sims = len(ALL_ARCHS) * len(CLIMATES)
    print(f"Initializing batch simulation: {total_sims} runs total.")
    print(f"Simulation step size set to: {SIMULATION_STEP} seconds.")
    
    tasks = []
    count = 0
    for arch_val_str in ALL_ARCHS:
        arch_value = float(arch_val_str)
        for climate in CLIMATES:
            count += 1
            tasks.append((arch_value, climate, count, total_sims))
            
    num_cores = multiprocessing.cpu_count()
    print(f"Running simulations in parallel using {num_cores} CPU cores...")
    
    # Track active running tasks and completion status in the main thread
    futures_map = {}
    completed_runs = 0
    
    with ProcessPoolExecutor(max_workers=num_cores) as executor:
        for t in tasks:
            arch_value, climate, count, total = t
            arch_sub_dir = f"Architecture_{int(arch_value)}"
            output_csv = os.path.join(RESULTS_BASE_DIR, arch_sub_dir, f"{climate}.csv")
            
            # Skip queueing if the file exists on disk already
            if os.path.exists(output_csv):
                print(f"[STATUS] Skipping (CSV already exists): Arch {arch_value} | Climate {climate} ({count}/{total})")
                completed_runs += 1
                continue
                
            print(f"[STATUS] Queueing: Arch {arch_value} | Climate {climate} ({count}/{total})")
            future = executor.submit(run_single_simulation, t)
            futures_map[future] = t
            
        if futures_map:
            for f in as_completed(futures_map):
                arch_value, climate, count, total = futures_map[f]
                completed_runs += 1
                try:
                    res = f.result()
                    print(f"[STATUS] ({completed_runs}/{total_sims}) Completed: {res}")
                except Exception as e:
                    print(f"[STATUS] ({completed_runs}/{total_sims}) Exception in Arch {arch_value} | Climate {climate}: {e}")
        else:
            print(f"[STATUS] All {total_sims} simulations already completed (CSVs exist on disk).")

# ----------------- KPI Data Extraction -----------------

def extract_data_worker(arch_id, folder_name):
    results = []
    for climate in CLIMATES:
        file_path = os.path.join(RESULTS_BASE_DIR, f"Architecture_{arch_id}", f"{climate}.csv")
        if not os.path.exists(file_path):
            file_path = os.path.join(RESULTS_BASE_DIR, folder_name, f"{climate}.csv")
        if not os.path.exists(file_path):
            file_path = file_path.replace(".csv", ".CSV")
            
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
                bins = ASHRAE_AIR_BINS if arch_id in AIR_ARCH_IDS else ASHRAE_LIQUID_BINS
                compliance = get_compliance(t_vals, bins)

                results.append({
                    'ArchID': int(arch_id), 'ArchName': ARCH_NAME_MAP.get(int(arch_id)),
                    'Climate': climate, 'PUE': pue, 'WUE': wue, 'CUE': cue,
                    'FC': fc_pct, 'PMC': pmc_pct, 'FMC': fmc_pct,
                    **compliance
                })
            except Exception as e: 
                print(f"Error reading/processing {file_path}: {e}")
    return results

# ----------------- Visualizations -----------------

def plot_comprehensive_kpis(df):
    plt.style.use('seaborn-v0_8-white')
    fig, axes = plt.subplots(3, 2, figsize=(10, 7), sharex=True)
    plt.subplots_adjust(hspace=0.08, wspace=0.18)
    
    df = df.sort_values('ArchID')
    all_arch_ids = range(1, 13)
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
    plt.savefig("KPI_Comprehensive_Analysis.svg", dpi=300, bbox_inches='tight')
    print("KPI_Comprehensive_Analysis charts saved successfully.")

def plot_kpi_circle_analysis(df):
    plt.style.use('seaborn-v0_8-white')
    metrics = ['PUE', 'WUE', 'CUE']
    summary_df = df.groupby('ArchID')[metrics].mean().reindex(range(1, 13))
    
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
        colors = plt.cm.get_cmap('tab20', 12)
    
    for i, arch_id in enumerate(range(1, 13)):
        vals = norm_df.loc[arch_id].values
        color = colors(i)
        
        # Calculate Centroid
        x_coords = vals * np.cos(angles)
        y_coords = vals * np.sin(angles)
        centroid_r = np.sqrt(np.mean(x_coords)**2 + np.mean(y_coords)**2)
        centroid_theta = np.arctan2(np.mean(y_coords), np.mean(x_coords))
        
        # Plot Triangle Outline
        plot_angles = np.append(angles, angles[0])
        plot_vals = np.append(vals, vals[0])
        ax.plot(plot_angles, plot_vals, color=color, linewidth=1.8, alpha=0.2, zorder=i+1, linestyle=":")

        # Plot Centroid
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
    plt.savefig("KPI_Efficiency_Circle_Triangles.svg", dpi=300, bbox_inches='tight')
    print("KPI_Efficiency_Circle_Triangles charts saved successfully.")

# ----------------- Main Coordinator -----------------

def main():
    # 1. Run batch simulations in parallel
    run_fmu_simulation()
    
    # 2. Extract results & Plot KPIs in parallel
    print("\nStarting parallel KPI data extraction and plotting...")
    if not os.path.exists(RESULTS_BASE_DIR):
        print(f"Results path not found: {RESULTS_BASE_DIR}")
        return
        
    folders = [f for f in os.listdir(RESULTS_BASE_DIR) if os.path.isdir(os.path.join(RESULTS_BASE_DIR, f))]
    all_map = { "".join(filter(str.isdigit, f)): f for f in folders if "".join(filter(str.isdigit, f)) in ALL_ARCHS }
    
    all_data = []
    num_cores = multiprocessing.cpu_count()
    print(f"Extracting metrics using {num_cores} CPU cores...")

    with ProcessPoolExecutor(max_workers=num_cores) as executor:
        futures = [executor.submit(extract_data_worker, aid, folder) for aid, folder in all_map.items()]
        for f in as_completed(futures):
            res = f.result()
            if res: 
                all_data.extend(res)
            
    if all_data:
        df_results = pd.DataFrame(all_data)
        plot_comprehensive_kpis(df_results)
        plot_kpi_circle_analysis(df_results)
        print("\nAll simulations and plotting completed successfully!")
    else:
        print("No valid simulation data found in results directory to plot.")

if __name__ == "__main__":
    multiprocessing.freeze_support()
    main()
