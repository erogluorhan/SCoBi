# SCoBi
[S]ignals of Opportunity [Co]herent [Bi]static Scattering Simulator


Repository contains:

1) source: Includes the following directories:
   
   a) input: The input contains:
   
      i) system: Simulation input files with the extension ".mat". For detailed information, please refer to the SCoBi user's manual.
      
      ii) configuration:  Configuration Input files (Excel ".xlsx"). For detailed information, please refer to the SCoBi user's manual.
      
      ii) Rx_antenna_pattern:  Receiver Antenna Pattern files (Excel ".xlsx"). For detailed information, please refer to the SCoBi user's manual.
      
      iii) vegetation:  Vegetation Input files (Excel ".xlsx"). For detailed information, please refer to the SCoBi user's manual.


      input directory contains two different sample inputs ( i.e. Paulownia for homogenous vegetation, and Corn for virtual vegetation)
      
   b) lib: It contains all the source codes of the SCoBi.

      To run SCoBi, basically:
      
      1) The function "runSCoBi.m" should be called.
      
      2) source/input should be ready.
	  
	  For further details about running simulations, please refer to SCoBi user's manual.
      
      IMPORTANT: When SCoBi is run, it generates the appropriate simulation outputs under "sims" directory, which MUST be IGNORED from repository commits. In other words, every developer or user of the repository should keep his/her "sims" directory LOCAL.


2) docs: Includes the following documents:
   
   a) manuals: The current SCoBi manuals are:
   
      i) User's Manual
      
      ii) Developer's Manual
      
      
3) design: Involves the architectural design of the SCoBi software that mostly consists of UML diagrams in Enterprise Architect.
