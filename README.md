# AmNat2022

README 

Harrison T, R Winfree, MA Genung. Price equations for understanding the response of ecosystem function to community change. The American Naturalist, in press.

Contact, responsible for data and code: Mark Genung, mark.genung@louisiana.edu

Summary: We present a novel, abundance-based version of the ecological Price equation in both discrete and continuous forms and explain the similarities and differences between this method and a related, previously developed richness-based method. We also present new empirical techniques for applying the Price equation to ecological data. The ecological Price equations derived here complement existing approaches, and together offer BEF researchers analytical tools and a unifying framework for studying biodiversity-ecosystem function in observational community data.
 
AmNat_Code_Feb22
Main R code for the paper. Produces all the results and figures described in the manuscript and appendices. Requires reading five .csv files, which are described below. The workflow for this code is simple, just run it from start to finish.

watermelon_com.csv
Site-species abundance matrix. Species are wild bees, sites are watermelon farms.

watermelon_com.csv
Site-species function matrix. Species are wild bees, sites are watermelon farms. Function is per-capita total grains of pollen deposited.

stream_inverts.csv
Site-taxa abundance matrix. Taxa are groups of stream invertebrate species, sites are streams. More details are available from the original data collectors at: https://datadryad.org/stash/dataset/doi:10.5061%2Fdryad.v6g985s

stream_invert_size.csv
Site-taxa function matrix. Taxa are groups of stream invertebrate species, sites are streams. Function is per-capita biomass. More details are available from the original data collectors at: https://datadryad.org/stash/dataset/doi:10.5061%2Fdryad.v6g985s

stream_data.csv
Site environmental characteristics for the stream invertebrate data. Site = steam where sampling occurred, A = abundance of invertebrates, s = species richness of invertebrates, cwm.z = community weighted mean function (biomass), f = total function (biomass), dist = disturbance, in this case anthropogenic water pollution

The analysis was run in R version 4.0.5 (2021-03-31) on the following platform: x86_64-apple-darwin17.0 (64-bit), running under: macOS Catalina 10.15.7. Required packages are included in the main code file.
![image](https://user-images.githubusercontent.com/8396183/153960351-3f0ac71c-7e46-4e4e-8fca-389798ab7069.png)

