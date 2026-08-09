# Software Documentation

**Author:** Viswanathan Ganesh  
**Date:** August 2026

---

## Abstract

This document provides comprehensive technical documentation for the **Data Center Thermal Management and Cooling Architecture Simulation Framework**. The framework couples high-fidelity dynamic thermal-fluid physics encapsulated within Functional Mock-up Units (FMUs) following the Functional Mock-up Interface (FMI 2.0) standard with multi-core Python execution pipelines. It models, simulates, and evaluates **12** distinct data center cooling topologies across **19** ASHRAE climate zones (0A through 8). Key performance metrics include Power Usage Effectiveness (PUE), Water Usage Effectiveness (WUE), Carbon Usage Effectiveness (CUE), operational mode distributions, and thermal compliance against ASHRAE environmental standards.

---

## Table of Contents

1. [System Overview & Core Capabilities](#1-system-overview--core-capabilities)
2. [Installation, Dependencies & Platform Compatibility](#2-installation-dependencies--platform-compatibility)
3. [Repository Directory Structure](#3-repository-directory-structure)
4. [Architecture & Operational Workflow](#4-architecture--operational-workflow)
5. [Execution Scripts: Full Batch vs. Benchmark Run](#5-execution-scripts-full-batch-vs-benchmark-run)
6. [Configuration Parameters & Sizing Criteria](#6-configuration-parameters--sizing-criteria)
7. [Input File Specifications & Data Schemas](#7-input-file-specifications--data-schemas)
8. [System Topology & Component Integration Matrix](#8-system-topology--component-integration-matrix)
9. [Mathematical & Physical Formulations](#9-mathematical--physical-formulations)
10. [Environmental Compliance Standards](#10-environmental-compliance-standards)
11. [Execution Workflow & Error Handling Safeguards](#11-execution-workflow--error-handling-safeguards)
12. [Output Artifacts & Publications](#12-output-artifacts--publications)

---

## 1. System Overview & Core Capabilities

The Data Center Thermal Management and Cooling Architecture Simulation Framework evaluates energy efficiency, water consumption, carbon footprint, and thermal compliance across **12** distinct data center cooling architectures operating under variable meteorological conditions defined by **19** ASHRAE climate zones (0A through 8).

The framework couples dynamic building physics models encapsulated inside a unified Functional Mock-up Unit (`DataCenterFMU.fmu`) via the Functional Mock-up Interface (FMI 2.0) co-simulation standard with high-performance Python computing pipelines. It orchestrates parallel parametric simulations, preprocesses meteorological boundary conditions, ingests time-varying IT load profiles, evaluates real-time thermodynamic metrics, and generates publication-quality key performance indicator (KPI) vector visualizations.

The repository provides two primary execution entry points:

1. **`master_run.py`** — Full multi-climate parametric batch execution engine configured across all 19 ASHRAE climate zones (0A to 8).
2. **`Example.py`** — A lightweight, targeted benchmark execution script configured specifically for ASHRAE Climate Zone 1A (Miami, Florida), serving as an illustrative validation run matching the primary reference scenario described in journal publications.

### 1.1 Key Features

- **Parallel Simulation Execution:** Dynamically distributes execution across multi-core CPU architectures using Python's `ProcessPoolExecutor` to accelerate multi-year parametric runs.
- **Single-Zone Validation Engine (`Example.py`):** Provides a rapid illustrative simulation run for Climate Zone 1A (Miami, FL) to reproduce journal article figures and benchmark results.
- **Dynamic IT Load Enforcement:** Ingests custom time-series IT workloads (`it_profile.csv`) and enforces strict physical isolation rules based on system topology.
- **C-Runtime Diagnostic Interception:** Captures lower-level solver callbacks and logs detailed failure diagnostics into `simulation_errors.log`.
- **OS-Level Stream Redirection:** Silences Newton-Raphson solver convergence warnings from compiled FMU C-DLLs during parallel execution.
- **Publication-Ready Visualizations:** Exports multi-panel KPI vectors, normalized radar charts, and equipment-level energy breakdown diagrams in vector graphics formats (`.pdf` and `.svg`).

---

## 2. Installation, Dependencies & Platform Compatibility

### 2.1 Prerequisites & Dependencies

The execution framework requires **Python 3.8 or higher** along with standard scientific computing packages:

```bash
pip install numpy pandas matplotlib fmpy
```

### 2.2 Platform Compatibility

The framework uses `FMPy` to instantiate and co-simulate FMUs. The primary model artifact (`DataCenterFMU.fmu`) contains pre-compiled dynamic-link libraries (C-DLLs) targeting **Win32** and **Win64** host operating systems.

---

## 3. Repository Directory Structure

```
.
├── Model/                          # Compiled FMU model directory
│   └── DataCenterFMU.fmu           # Unified cooling loop FMU model
├── Climate Files/                  # Weather profile files (.mos format)
│   ├── 0A.mos
│   ├── 1A.mos                      # Primary benchmark profile (Miami, FL)
│   ├── ...
│   └── 8.mos
├── it_profile.csv                  # Time-varying IT load profile
├── master_run.py                   # Multi-climate batch execution script
├── Example.py                      # Illustrative benchmark script (Zone 1A)
├── simulation_errors.log           # Diagnostic log capturing C-runtime warnings
├── Documentation.pdf               # Repository overview
└── results/                        # Execution CSV outputs directory
    ├── Architecture_1/
    │   ├── 1A.csv
    │   └── ...
    └── Architecture_12/
```

---

## 4. Architecture & Operational Workflow

The simulation execution pipeline follows a **five-stage modular architecture**:

1. **Meteorological & Load Preprocessing:** Parses EnergyPlus/Modelica `.mos` weather files, derives wet-bulb temperatures using Stull's empirical relation, ingests dynamic IT load profiles from `it_profile.csv`, and performs temporal linear interpolation to the communication step size.

2. **Parallel FMU Orchestration Engine:** Distributes tasks via `ProcessPoolExecutor`, injects dynamic setpoint vector inputs, enforces relative solver tolerances (10^-4) to prevent hydronic flow chatter, and redirects OS-level stdio/stderr streams to suppress C-DLL solver warnings.

3. **C-Runtime Exception & Error Interceptor:** Intercepts internal FMU C-logger callbacks (status >= 1), categorizes numerical non-convergence, temperature boundary exceedances, or mass balance warnings, and appends diagnostic traces to `simulation_errors.log`.

4. **Thermodynamic & KPI Evaluation Engine:** Calculates pointwise Water Usage Effectiveness (WUE), Power Usage Effectiveness (PUE), Carbon Usage Effectiveness (CUE), classifies operating modes (Free Cooling, Partial Free Cooling, Full Mechanical Cooling), and performs ASHRAE compliance binning.

5. **Visualization & Report Generation:** Generates multi-panel KPI summaries, normalized efficiency radar plots, and equipment-level energy breakdown bar charts exported as publication-ready vector PDFs.

---

## 5. Execution Scripts: Full Batch vs. Benchmark Run

The repository contains two execution scripts designed for distinct analytical contexts:

| Parameter / Feature | `master_run.py` | `Example.py` |
|---|---|---|
| **Target Purpose** | Multi-climate parametric analysis | Single-zone benchmark validation |
| **Climate Zones Simulated** | 19 Zones (0A to 8) | 1 Zone (1A — Miami, Florida) |
| **Primary Context** | Regional thermodynamic sweep | Journal article case study |
| **Total Simulation Runs** | 228 (12 Archs x 19 Climates) | 12 (12 Archs x 1 Climate) |
| **Output Artifacts** | `results/Architecture_X/*.csv` + 3 KPI PDF Graphics | `results/Architecture_X/1A.csv` + 3 KPI PDF Graphics |

---

## 6. Configuration Parameters & Sizing Criteria

### 6.1 IT Infrastructure Sizing & Architectural Isolation Rules

The cooling equipment in the unified FMU model is sized for a maximum nominal IT load **P_IT,max = 500 kW** (500,000 W).

To model distinct system boundaries without structural re-compilation, the orchestrator applies strict load isolation rules:

- **Air-Cooled Systems (Types 1-6):** The framework sets the liquid cooling loop load to zero (P_IT,Liq = 0.0 W), assigning the full load profile to the air loop (P_IT,Air = P_IT).
- **Liquid-Cooled Systems (Types 7-12):** The framework isolates the air loop by setting P_IT,Air = 0.0 W, routing the full IT load to the liquid loop (P_IT,Liq = P_IT).

### 6.2 Simulation Timestep & Setpoint Parameters

| Parameter | Value |
|---|---|
| **Simulation Communication Step (delta_t_sim)** | 3600 s (1 hour) |
| **Output Sampling Interval (delta_t_out)** | 3600 s |
| **Simulation Horizon** | 365 days (8760 hours) |
| **Chilled Water Supply Setpoint (T_CHW,sup)** | 15.0 degC (288.15 K) |
| **CDU Supply Setpoint (T_CDU,sup)** | 27.0 degC (300.15 K) |
| **Air Supply Setpoint (T_air,sup)** | 18.0 degC (291.15 K) |
| **Air Return Setpoint (T_air,ret)** | 25.0 degC (298.15 K) |
| **Liquid Return Setpoint (T_liq,ret)** | 32.0 degC (305.15 K) |
| **Freeze Protection Threshold (T_freeze)** | 5.0 degC (278.15 K) |
| **Approach Temperature (delta_T_app)** | 6.0 K (Dry Systems) / Variable (Wet Systems) |

---

## 7. Input File Specifications & Data Schemas

### 7.1 Meteorological Profile Requirements (`.mos`)

Weather profile files reside in `Climate Files/` and follow the Modelica tab-separated `.mos` structure. Meteorological profiles can be retrieved from the [EnergyPlus Weather Database](https://energyplus.net/weather) and formatted into four tab-delimited columns:

1. **Column 1:** Annual elapsed time in seconds (t in [0, 31,536,000]).
2. **Column 2:** Ambient dry-bulb temperature (T_db, in degC).
3. **Column 3:** Relative humidity (RH, in %, 0-100).
4. **Column 4:** Ambient wet-bulb temperature (T_wb, in degC) or dynamically derived psychrometrically.

### 7.2 Dynamic IT Load Profile Schema (`it_profile.csv`)

Custom dynamic IT loads are provided via `it_profile.csv` located at the project root:

| Column Name | Data Type | Units | Description |
|---|---|---|---|
| `time` | Float | Seconds (s) | Annual timestamp (0 to 31,536,000 s) |
| `ITAir` | Float | Watts (W) | Heat load for air-cooled loop (<= 500,000 W) |
| `ITLiq` | Float | Watts (W) | Heat load for liquid-cooled loop (<= 500,000 W) |

---

## 8. System Topology & Component Integration Matrix

The 12 evaluated topologies are divided into **Air-Cooled Systems (Types 1-6)** and **Liquid-Cooled Systems (Types 7-12)**:

| Arch ID | Medium | Chiller | WSE | WCT | DCT | CW Pump | CHW | AHU | CDU |
|---|---|---|---|---|---|---|---|---|---|
| **Type 1**  | Air    | --  | X | X | --  | X | X | X | --  |
| **Type 2**  | Air    | X  | X | X | --  | X | X | X | --  |
| **Type 3**  | Air    | X  | --  | X | --  | X | X | X | --  |
| **Type 4**  | Air    | --  | X | --  | X | X | X | X | --  |
| **Type 5**  | Air    | X  | X | --  | X | X | X | X | --  |
| **Type 6**  | Air    | X  | --  | --  | X | X | X | X | --  |
| **Type 7**  | Liquid | --  | X | X | --  | X | X | --  | X |
| **Type 8**  | Liquid | X  | X | X | --  | X | X | --  | X |
| **Type 9**  | Liquid | X  | --  | X | --  | X | X | --  | X |
| **Type 10** | Liquid | --  | X | --  | X | X | X | --  | X |
| **Type 11** | Liquid | X  | X | --  | X | X | X | --  | X |
| **Type 12** | Liquid | X  | --  | --  | X | X | X | --  | X |

> **Legend:** WSE = Water-Side Economizer, WCT = Wet Cooling Tower, DCT = Dry Cooling Tower / Dry Cooler, CW Pump = Condenser Water Pump, CHW = Chilled Water Loop, AHU = Air Handling Unit, CDU = Coolant Distribution Unit.

---

## 9. Mathematical & Physical Formulations

### 9.1 Psychrometric Calculations

When wet-bulb temperature $T_{\text{wb}}$ is not directly supplied by the meteorological profile, it is calculated using **Stull's empirical psychrometric formulation**:

$$T_{\text{wb}} = T_{\text{db}} \cdot \arctan\!\left(0.151977\sqrt{RH + 8.313659}\right) + \arctan(T_{\text{db}} + RH) - \arctan(RH - 1.676331) + 0.00391838 \cdot RH^{1.5} \cdot \arctan(0.023101 \cdot RH) - 4.686035$$

Saturation vapor pressure $p_{\text{sat}}(T)$ (kPa) is derived via **Tetens' relation**:

$$p_{\text{sat}}(T) = 0.61121 \cdot \exp\!\left(\frac{(18.678 - T/234.5)\cdot T}{257.14 + T}\right)$$

For standard atmospheric pressure $P_{\text{atm}} = 101.325\text{ kPa}$, **humidity ratio** $w$ (kg_water / kg_dry air) is:

$$w = \frac{h_{\text{fg,wb}} \cdot w_{\text{s,wb}} - 1.006\,(T_{\text{db}} - T_{\text{wb}})}{h_{\text{fg,wb}} + 1.86\,T_{\text{db}} - 4.186\,T_{\text{wb}}}$$

where:

$$w_{\text{s,wb}} = 0.622 \cdot \frac{p_{\text{sat}}(T_{\text{wb}})}{P_{\text{atm}} - p_{\text{sat}}(T_{\text{wb}})}$$

$$h_{\text{fg,wb}} = 2501 - 2.326 \cdot T_{\text{wb}} \quad [\text{kJ/kg}]$$

**Specific enthalpy** $h$ (kJ/kg):

$$h = 1.006\,T_{\text{db}} + w\,(2501 + 1.86\,T_{\text{db}})$$

---

### 9.2 Dynamic Water Usage Effectiveness (WUE) Engine

Makeup water consumption accounts for evaporation, blowdown, and drift:

**1. Heat Rejection Rate:**

$$\dot{Q}_{\text{rej}} = \dot{m}_{\text{CW}} \cdot c_{p,\text{water}} \cdot (T_{\text{CW,ret}} - T_{\text{CW,sup}})$$

where $c_{p,\text{water}} = 4.186\text{ kJ/(kg·K)}$.

**2. Evaporation Rate** — Air exiting the cooling tower is assumed saturated at $T_{\text{exit}} = T_{\text{wb}} + 3.0\,^\circ\text{C}$:

$$\dot{m}_{\text{evap}} = \dot{Q}_{\text{rej}} \cdot \left(\frac{\Delta w}{\Delta h}\right)$$

**3. Total Makeup Water Flow Rate:**

$$\dot{m}_{\text{makeup}} = \dot{m}_{\text{evap}} + \frac{\dot{m}_{\text{evap}}}{COC - 1} + (\dot{m}_{\text{CW}} \cdot \delta_{\text{drift}})$$

where $COC = 5.0$ (Cycles of Concentration) and $\delta_{\text{drift}} = 0.0002$ (0.02% drift loss rate).

**4. Pointwise WUE (L/kWh):**

$$\text{WUE}(t) = \frac{\dot{m}_{\text{makeup}}(t) \cdot 3600}{P_{\text{IT}}(t)}$$

---

### 9.3 PUE, CUE, and Energy Consumption Metrics

Total HVAC power consumption:

$$P_{\text{HVAC}} = P_{\text{Chiller}} + P_{\text{CT,Fan}} + P_{\text{Pump,CW}} + P_{\text{Pump,CHW}} + P_{\text{AHU,Fan}} + P_{\text{Pump,CDU}} + P_{\text{DCT,Fan}}$$

**Power Usage Effectiveness (PUE)** and **Carbon Usage Effectiveness (CUE)**:

$$\text{PUE}(t) = 1.0 + \frac{P_{\text{HVAC}}(t)}{P_{\text{IT}}(t)}$$

$$\text{CUE}(t) = \text{PUE}(t) \times \text{CI}_{\text{climate}}$$

where $\text{CI}_{\text{climate}}$ represents grid carbon intensity (kg CO2/kWh).

---

### 9.4 Operational Mode Classification Logic

Operating modes are tracked dynamically via Chiller valve modulation ($y_{\text{Chi}} \in [0, 1]$) and WSE valve modulation ($y_{\text{WSE}} \in [0, 1]$):

$$\text{Mode}(t) = \begin{cases} \text{Full Free Cooling (FC)}, & \text{if } y_{\text{Chi}} = 0 \text{ and } y_{\text{WSE}} = 1 \\ \text{Partial Free Cooling (PMC)}, & \text{if } y_{\text{Chi}} = 1 \text{ and } y_{\text{WSE}} = 1 \\ \text{Full Mechanical Cooling (FMC)}, & \text{if } y_{\text{Chi}} = 1 \text{ and } y_{\text{WSE}} = 0 \end{cases}$$

---

## 10. Environmental Compliance Standards

### 10.1 ASHRAE Guidelines for Air-Cooled Systems

Air compliance evaluates supply air temperature $T_{\text{air,sup}}$ (degC):

| Class | Temperature Range |
|---|---|
| **Recommended** | 18.0 degC <= T_air,sup <= 27.0 degC |
| **Class A1** | 15.0 degC <= T_air,sup <= 32.0 degC |
| **Class A2** | 10.0 degC <= T_air,sup <= 35.0 degC |
| **Class A3** | 5.0 degC <= T_air,sup <= 40.0 degC |
| **Class A4** | 5.0 degC <= T_air,sup <= 45.0 degC |

### 10.2 ASHRAE Guidelines for Liquid-Cooled Systems

Liquid compliance evaluates CDU supply coolant temperature $T_{\text{CDU,sup}}$ (degC):

| Class | Temperature Range |
|---|---|
| **Class W17** | 2.0 degC <= T_CDU,sup <= 17.0 degC |
| **Class W27** | 2.0 degC <= T_CDU,sup <= 27.0 degC |
| **Class W32** | 2.0 degC <= T_CDU,sup <= 32.0 degC |
| **Class W40** | 2.0 degC <= T_CDU,sup <= 40.0 degC |
| **Class W45** | 2.0 degC <= T_CDU,sup <= 45.0 degC |

**Compliance percentage** $C_{\text{bin}}$ over total timesteps $N$:

$$C_{\text{bin}} = \left(\frac{1}{N} \sum_{t=1}^{N} \mathbf{1}\!\left(T_{\text{low}} \le T(t) \le T_{\text{high}}\right)\right) \times 100\%$$

---

## 11. Execution Workflow & Error Handling Safeguards

### 11.1 Command-Line Execution Commands

**Option A — Illustrative Benchmark Run (Zone 1A - Miami, FL):**

```bash
python Example.py
```

**Option B — Complete Multi-Climate Batch Simulation:**

```bash
python master_run.py
```

### 11.2 Robustness & Diagnostics

1. **OS-Level Stream Redirection:** Direct system calls (`stdout`=1, `stderr`=2) are suppressed at the OS layer during FMU solver execution to prevent buffer lockups during multi-core processing.
2. **Numeric Integration Tolerances:** Relative solver integration tolerance is set to 10^-4, suppressing high-frequency chattering and flow reversal instabilities in hydronic loops.
3. **C-Runtime Exception Interception:** Internal FMU C-logger callbacks capture status flags >= 1 and record explicit failure details in `simulation_errors.log`.

---

## 12. Output Artifacts & Publications

The execution scripts output detailed time-series CSV files under `results/Architecture_<ID>/<Climate>.csv` and automatically export high-resolution vector figures:

1. **`KPI_Comprehensive_Analysis.pdf`** — Multi-panel chart detailing PUE, WUE, CUE, operational mode distributions, and ASHRAE thermal compliance percentages.
2. **`KPI_Efficiency_Circle_Triangles.pdf`** — Polar-coordinate radar plot mapping normalized performance centroids for PUE, WUE, and CUE.
3. **`KPI_Equipment_Energy_Consumption.pdf`** — Stacked bar chart breaking down annual electrical energy consumption (MWh) across Chillers, Tower Fans, Dry Cooler Fans, Condenser Water Pumps, Chilled Water Pumps, AHU Fans, and CDU Pumps.
