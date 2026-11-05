R reproduction of the supplied one-carbon/Fe-S/nucleotide pathway schematic.

Run: Rscript draw_pathway.R /path/to/output
Requires base R with grid. Tested with R 4.3.3. svglite is optional and preserves text elements in the SVG.

Files: draw_pathway.R, vector PDF and SVG, 300 dpi PNG, label/edge coordinate tables, and R session information.
Edit label() and edge() definitions in the R script to change content or layout. Coordinates use a 1120 x 680 canvas with origin at the top left.

The biological labels and connections reproduce the user-supplied diagram; this deliverable is a drawing reconstruction rather than a new pathway validation.
