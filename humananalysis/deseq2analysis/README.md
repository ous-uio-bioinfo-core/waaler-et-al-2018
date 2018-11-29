# DESeq2  analysis

For the comparison of differences in effects between a group of 4 cell lines and another group of 9 cell lines an extra test was performed using DESeq2. The same test was performed using limma, however the design was somewhat complicated, and we were not sure that we got it right so this alternative test was done using DESeq2 where is was possible to model cell line as a fixed effect in the design (whereas it was set as a random effect in limma).
<br/>
<br/>
DESeq2 identifies several genes as differentially expressed, i.e the effect of the drug differ between the groups. Limma does not find these genes as DEG.
<br/>
<br/>
This scripts may be reproducible with some minor alterations.
<br/>
<br/>



