# SIM-BE
Fast multistage carcinogenesis simulation of neoplastic progression in Barrett's Esophagus 

[see accompanying R-Notebook (presently unfinished)]

This repository provides code for a superfast simulation of the clonal evolution in Barrett's esophagus. It simulates pairs of premaligant (Y) clones and associated  'first-passage-time' malignant (M) clones and records clone sizes and the origins of each Y-clone and (first) M-clone. The simulation starts with a determination of a random BE onset time and a random set of initiated progenitor cells in a BE segment for a given individual. Each initiated cell develops into a stochastic Y-clone which may or may not make a malignant transformation in the future. This is determined upfront for all initiated cells. Thus, there are two types of Y-clones: Y-clones that make a 'first-passage-time' (FPT) mutation (recorded in DF $Malign$), and those that don't (recorded on DF $Benign$). The latter will go to extinction with certainty and will never make a malignant transformation. For all pratical purposes, they can be neglected when the extinction probability is high. For speed, the two 'propagator' functions, $GENbenign$ and $GENmalign$, are set up to be called recursively, generating random clone sizes at 'screening' times $t_i$, with $t_i < t_{i+1}$.

Simulated observations are stored in two dataframes ($Benign$ and $Malign$) and updated in bulk as time progresses from BE onset to MaxAge:

for example, for $Malign$, by column \
1) maximum age since onset of BE (BEage) \
2) start: vector of births of initiated cells destined to develop into clones that undergo a malignant transformation (s.initiated.malignant) \
3) s: running time (counting from BE onset) \ 
4) Y: initiated clone sizes (Y-clones) \
5) Z: malignant clone size (M-clones) \
6) fpT: first passage times \
7) LGD, HGD, IMC, EAC logical indicators \

The simulation relies on just four R modules/functions. They use (global) parameters defined in *BE_MSCE1_modelParameters.R* and *BE_MSCE1_shared.R*.

**TC1**: simulation of first passage times (eg first malignant transformations) among $n$ evolving Y-clones using the 'one-stage' tumor survival function $S_1$ with $(Y(0)=1)$, applied indpendently to the $n$ progenitor cells. First passage times can be $\infty$ if a clone goes extinct before giving rise to a maligant transformation. 

others to be provided in R Notebook
