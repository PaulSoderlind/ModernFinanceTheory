# Introduction

This repository contains Julia code for a MSc course in Financial Theory at UNISG. 


# Instructions

1. Most files are jupyter notebooks. Click one of them to see it online. If GitHub fails to render the notebook or messes up the LaTeX in the Markdown cells, then use [nbviewer](https://nbviewer.jupyter.org/). Instructions: try to open the notebook at GitHub, copy the link and paste it in the address field of nbviewer.

2. Two pdf files contain print-outs of all notebooks. 

3. To download this repository, use the Download (as zip) in the Github menu. Otherwise, clone it.


# On the Files

1. ChapterNumber_Topic.ipynb are notebooks organised around different topics.

2. NotebooksAsPDFxxxx.pdf is a print-out of all notebooks. 

3. The folder Data contains data sets used in the notebooks, while the folder src contains .jl files with some functions also used in the notebooks.

4. The plots are in png format. If you want sharper plots, change `default(fmt = :png)` to `default(fmt = :svg)` in one of the top cells.

5. The current version is tested on Julia 1.12.

