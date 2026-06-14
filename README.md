# Data Center Cooling Model Simulation (FMI / Co-Simulation)

This repository contains a self-contained batch simulation and post-processing pipeline for data center cooling loop configurations. It utilizes Functional Mock-up Units (FMUs) via FMPy and evaluates cooling performance across multiple structural architectures.

> [!IMPORTANT]
> - **Equipment Sizing**: The cooling equipment in the provided FMU model is sized for a **maximum IT load of 1 MW** (1,000,000 W).
> - **Climate Configurations**: The pipeline is currently configured to run for **one climate file** (`DataCenter_0A`). Support for additional climate configurations will be explored in future updates.


## Key Features

- **Parallel Simulation Execution**: Distributes the simulation runs dynamically across all available CPU cores using a multiprocessing pool to drastically reduce total runtime.
- **Dynamic IT Load Profiles**: Accepts a custom, time-varying IT load profile (`it_profile.csv`) and automatically interpolates values to match simulation step intervals.
- **Architectural Rules & Validation**:
  - **Type 1 to 6 (Air-Cooled Only)**: Automatically forces liquid-cooling loop loads to zero (`ITLiq = 0.0 W`).
  - **Type 7 to 12 (Liquid-Cooled Only)**: Automatically forces air-cooling loop loads to zero (`ITAir = 0.0 W`).
- **OS-Level Solver Warning Redirection**: Silences Newton solver convergence warnings from the compiled C-runtime DLLs for clean console logging.
- **Interactive KPI Visualizations**: Exports publication-quality vector charts (`.pdf` and `.svg`):
  - **KPI Comprehensive Analysis**: Compares PUE, WUE, CUE, operating mode states (Free Cooling, Partially Mechanically Cooled, Fully Mechanically Cooled), and ASHRAE compliance bins.
  - **KPI Efficiency Circle Triangles**: Normalizes and maps PUE, WUE, and CUE on a polar axis.


---

## Cooling System Architecture Types

The batch simulation pipeline supports the following 12 cooling system architecture configurations (Types 1–6 are Air-cooled; Types 7–12 are Liquid-cooled):

| Cooling Medium | Cooling System | Chiller | WSE | WCT | DCT | CW Pump | CHW | AHU | CDU |
| :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **Air** | **Type 1** | | X | X | | X | X | X | |
| | **Type 2** | X | X | X | | X | X | X | |
| | **Type 3** | X | | X | | X | X | X | |
| | **Type 4** | | X | | X | X | X | X | |
| | **Type 5** | X | X | | X | X | X | X | |
| | **Type 6** | X | | | X | X | X | X | |
| **Liquid** | **Type 7** | | X | X | | X | X | | X |
| | **Type 8** | X | X | X | | X | X | | X |
| | **Type 9** | X | | X | | X | X | | X |
| | **Type 10**| | X | | X | X | X | | X |
| | **Type 11**| X | X | | X | X | X | | X |
| | **Type 12**| X | | | X | X | X | | X |

*Legend: WSE = Water Side Economizer, WCT = Water Cooling Tower, DCT = Dry Cooling Tower (Dry Cooler), CW Pump = Condenser Water Pump, CHW = Chilled Water Loop, AHU = Air Handling Unit, CDU = Coolant Distribution Unit.*

---

## Repository Structure

```
├── DataCenter_0A.fmu        # Compiled FMU model configuration
├── it_profile.csv          # Custom IT load profiles input (time, ITAir, ITLiq)
├── master_run.py           # Self-contained simulation and plotting script
├── README.md               # Repository documentation
└── results/                # Output directory containing generated CSVs (created automatically)
```

---

## Installation & Requirements

Ensure you have Python 3.8+ installed along with the required libraries:

```bash
pip install numpy pandas matplotlib fmpy
```

*Note: FMPy requires a C-compiler or pre-compiled binaries to execute FMU Co-Simulations. The FMUs in this repository contain pre-compiled DLLs for Win32 and Win64 platforms.*

---

## Usage

### 1. Configure Setpoints (Optional)
Open **[master_run.py](file:///d:/SoftwareX/master_run.py)** and modify the global temperature setpoints in Celsius at the top:
```python
# Temperature Setpoints in Celsius
CHW_SUP_SET_C = 15.0  # Chilled Water Supply Setpoint (°C)
CDU_SUP_SET_C = 22.0  # CDU Supply Setpoint (°C)
```

### 2. Define Custom IT Profiles (Optional)
Modify **[it_profile.csv](file:///d:/SoftwareX/it_profile.csv)** to input your own IT load profiles. It takes:
- `time`: Time elapsed from the beginning of the year (in seconds).
- `ITAir`: Heat load in air-cooled server rooms (in Watts).
- `ITLiq`: Heat load in liquid-cooled systems (in Watts).

### 3. Run the Pipeline
Execute the master script:
```bash
python master_run.py
```

### 4. Output Results
Upon execution, the script will:
1. Create a `results/` folder containing subfolders for each architecture type (e.g. `results/Architecture_1/DataCenter_0A.csv`).
2. Generate and save two visualization plots:
   - **`KPI_Comprehensive_Analysis.pdf` / `.svg`**
   - **`KPI_Efficiency_Circle_Triangles.pdf` / `.svg`**
