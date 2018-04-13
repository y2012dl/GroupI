1. Execution
PyMUS v2.0 is optimized for Windows 7 and 8.1. The Main Window of PyMUS v2.0 can be launched in two ways. If you like to use the source code, please run "GUI.py". Otherwise, you can run the executable file of "PyMUS.exe" depending on the version of your operating system (32-bit or 64-bit).

2. Dependencies
- Python (ver 2.7)
- Pandas (ver 0.15.2)
- Numpy (ver 1.9.2)
- Scipy (ver 0.15.1)
- Matplotlib (ver 1.4.3)
- PyQt4 (ver 4.10.4)

3. Software Usage
PyMUS v2.0 allows the simulations to be fully operated using a graphical user interface (GUI). The GUI was designed to allow for a generic computational procedure to be performed using modeling and simulation approaches. The Main Window of PyMUS v2.0 consists of one state window and six buttons for controlling simulation of motor unit system. A typical procedure for use of the PyMUS v2.0 is as follows:

[STEP 1] Model Selection
 - Using ¡°Please Select!¡± drop-down menu, the target model to be simulated is chosen among the motoneuron, muscle unit and motor unit models.
 - The message on the result of STEP1 appears on the state window.

[STEP 2] Parameter Setting
 - Using ¡°Model Parameter Settings¡± button, the Model Parameters window is popped up to provide the GUI interfaces for setting the model parameter values manually or automatically by importing pre-determined data into the software.
 - When ¡°OK¡± or ¡°Apply¡± button is pushed, STEP2 is completed along with the result message on the state window. 

[STEP 3] Simulation Setting
  - Using ¡°Simulation Condition Settings¡± button, the Simulation Conditions window is generated to provide the GUI interfaces for setting the simulation time, display quality (sample time and plotting interval) and initial values of the model equations. 
  - When ¡°OK¡± or ¡°Apply¡± button is pushed, STEP3 is completed along with the result message on the state window.

[STEP 4] Input Setting
  - Using ¡°Input Signal Settings¡± button, the Input Signals window is produced to allow the type and protocol for injecting the input signals into the model to be selected.
  - When ¡°Generate¡± button is pushed, input signal data is produced and displayed on a separate figure window along with the result message on the state window. 
  - When ¡°OK¡± or ¡°Apply¡± button is pushed, STEP4 is completed along with the result message on the state window.

[STEP 5] Output Setting
- Using ¡°Output Signal Settings¡± button, the Output Signals window is popped up to enable the output variables to be displayed and plotting options to be selected.
- To efficiently compare the multiple output variables simultaneously, PyMUS v2.0 allows for the output variables to be selected in the Output Signals window and displays these variables either individually on separate panels or together on the same panel.
 - When ¡°OK¡± or ¡°Apply¡± button is pushed, STEP5 is completed along with the result message on the state window. 

[STEP6] Run control
  - Using ¡°Run¡± button, the simulation can be started or stopped.
  - The messages on the result of STEP6 are displayed on the state window.

STEP2-STEP6 can mutually transition among each other maintaining all changes made in each step unless the target model is changed through STEP1. To support offline analyses, the simulation data including simulation time and all output variables in the Output Signals window are also saved in a separate file.