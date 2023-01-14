These data were used in the creation of the article "Tectonic advection of contacts enhances landscape transience," which is currently in review with Earth Surface Processes and Landforms.

Each set of simulations has an Excel input table. Not all of the simulations listed in the tables
were used in the article, however. Additionally, the set numbers do not match those in Table 1 
of the article.

Sets 1 and 2 use topographic advection rates of v = 0.5 mm/yr.
Sets 3 and 4 use topographic advection rates of v = 0.25 mm/yr.
Sets 5 and 6 use topographic advection rates of v = 0 mm/yr. Simulations in Sets 5 and 6 have different 
contact advection rates (vc): 0.5, 0.25, and 0 mm/yr.

Note that topographic advection rates (v) and contact advection rates (vc) are in the positive y 
direction (i.e., to the north). No advection accurs in the x direction (West to East).

For simulation Sets 1 through 4 (which include topographic advection), the simulations used include 
9, 10, and 11 (Kw/Ks values of 2, 5, and 10 and a contact dip of 90 degrees) as well as 21, 22, and 
23 (Kw/Ks values of 2, 5, and 10 and a contact dip of 30 degrees).

For simulation Sets 5 and 6 (which do not include topographic advection), th simulations used include:
9, 10, and 11 (Kw/Ks values of 2, 5, and 10 and a contact dip of 90 degrees with contact advection rate 
vc = 0.5 mm/yr); 13, 14, and 15 (Kw/Ks values of 2, 5, and 10 and a contact dip of 90 degrees with contact
advection rate vc = 0.25 mm/yr); 21, 22, and 23 (Kw/Ks values of 2, 5, and 10 and a contact dip of 30 degrees 
with contact advection rate vc = 0.5 mm/yr); 25, 26, and 27 (Kw/Ks values of 2, 5, and 10 and a contact dip
of 30 degrees with contact advection rate vc = 0.25 mm/yr); and 33, 34, and 35 (Kw/Ks values of 2, 5, and 10 
and a contact dip of 30 degrees with contact advection rate vc = 0 mm/yr).

Due to the large sizes of the output files, I have only included the last timestep for the initialization of
each simulation set. These initialized landscapes and the scripts included can be used to rerun each scenario.

In Sets 1-4,  the scenarios with contact dips of 90 degrees (9-11) can be run with the code 
"Run_scenarios_unlimited_weak_zone.m."
In Sets 1-4,  the scenarios with contact dips of 30 degrees (21-23) can be run with the code 
"Run_scenarios_unlimited_weak_zone_Dip.m."

In Sets 5-6,  the scenarios with contact dips of 90 degrees (9-11 and 13-15) can be run with the 
code "Run_scenarios_unlimited_weak_zone_no_advec.m."
In Sets 5-6,  the scenarios with contact dips of 30 degrees (21-23, 25-27, and 33-35) can be run 
with the code "Run_scenarios_unlimited_weak_zone_no_advec_Dip.m."

Just update the Scenario number within the codes mentioned above and ensure that the data for the initialized
landscape (in a .mat file) is present within the Output folder (in the same directory as the code and Excel table).
Within the Output folder, the .mat file should be in a folder called "Set#_Initialization," where # is replaced by 
the set number. These scenarios can take a very long time to run.

The edited TTLEM codes (which are required for these simulations) are included in the folder "Edited_TTLEM_Codes."
