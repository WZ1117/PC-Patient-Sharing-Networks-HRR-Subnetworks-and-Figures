# PC Patient-Sharing Networks (HRR Subnetworks & Demo)
This repo is for anyone who wants to try building and visualizing networks from patient–provider data. Includes R code to create provider–provider networks, color by specialty, and plot HRR subnetworks. You can use your own dataset or generate synthetic demo data.

## What’s here
- **make_network_and_plot.Rmd** – end-to-end: load → build bipartite/unipartite → attach attributes → plot HRR subnetworks
- **network_figures/** – output images (created on run)

## Expected columns
`patient_id, npi, provider_hrr, provider_specialty, pc_specialist_flag, (optional: patient_hrr, cancer_site, pc_flag, encounter_year)`

## Run
Open `make_network_and_plot.Rmd` in RStudio and Knit, or run the chunks sequentially.  
Outputs are saved to `network_figures/`.

## Notes
- Edges are weighted by **shared patients** between providers (from bipartite projection).
- Node color = specialty group; PC providers are framed in black and plotted as squares; non-PC providers are plotted as circles.
- We showe here a **schematic demo** figure rather than raw clinician-level plots for reference how the figures will look like.
