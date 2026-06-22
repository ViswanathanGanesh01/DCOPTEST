# Data Center Cooling Model Simulation (FMI / Co-Simulation)

This repository contains a self-contained batch simulation and post-processing pipeline for data center cooling loop configurations. It utilizes Functional Mock-up Units (FMUs) via FMPy and evaluates cooling performance across multiple structural architectures.

> [!IMPORTANT]
> - **Equipment Sizing**: The cooling equipment in the provided FMU model is sized for a **maximum IT load of 500 kW** (500,000 W).
> - **FMU Model Configuration**: The pipeline runs using a single, unified FMU model **`DataCenterFMU.fmu`** placed in the `Model/` directory, driven dynamically by external climate `.mos` files to evaluate performance.


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
├── Model/                  # Directory containing compiled FMU model configuration
│   └── DataCenterFMU.fmu   # Self-contained cooling loop FMU model
├── it_profile.csv          # Custom IT load profiles input (time, ITAir, ITLiq)
├── master_run.py           # Self-contained simulation and plotting script
├── README.md               # Repository documentation
├── results/                # Output directory containing generated CSVs (created automatically)
```

## Input File Requirements

### 1. Climate / Weather Files (`.mos` Format)
The simulation requires hourly weather profile files placed in the `Climate Files/` directory.
- **Source**: Weather files can be downloaded from the [EnergyPlus Weather Website](https://energyplus.net/weather).
- **Format**: The weather files must be in the standard Modelica weather format (`.mos` tab-separated files) containing:
  - Column 1: Time in seconds (from 0 to 31,536,000 for a full year).
  - Column 2: Dry bulb temperature in Celsius.
  - Column 3: Relative humidity (0–100%).
  - Column 4: Wet bulb temperature in Celsius (or calculated dynamically).
- **Conversion Note**: If you download a weather profile in another format (like `.epw`), it must be converted or reformatted into the four-column tab-separated `.mos` layout before execution.

### 2. IT Load Profile (`it_profile.csv`)
You can define custom time-varying IT loads using the `it_profile.csv` file in the root folder:
- **Columns**:
  - `time`: Time step in seconds from the start of the year (0 to 31,536,000).
  - `ITAir`: Heat load on the air-cooling loop (in Watts). Sized for up to 500,000 W (500 kW).
  - `ITLiq`: Heat load on the liquid-cooling loop (in Watts). Sized for up to 500,000 W (500 kW).
- **Note**: The pipeline automatically forces `ITAir` to zero for liquid architectures (Types 7–12), and `ITLiq` to zero for air architectures (Types 1–6).

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
1. Create a `results/` folder containing subfolders for each architecture type (e.g. `results/Architecture_1/1A.csv`).
2. Generate and save the visualization plots:
   - **`KPI_Comprehensive_Analysis.pdf`**
   - **`KPI_Efficiency_Circle_Triangles.pdf`**
   - **`KPI_Equipment_Energy_Consumption.pdf`**
3. **Error Logging**: If any cooling topology fails during simulation, the script dynamically intercepts the C-level FMU warnings/errors, categorizes the failure reason (e.g., *Temperature above/below limit*, *Mass fraction bounds*, or *Numerical stiffness*), and appends it to **`simulation_errors.log`** in the root directory.
