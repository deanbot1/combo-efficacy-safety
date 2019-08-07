# License & Disclaimer
MIT License

Copyright (c) 2019 Millennium Pharmaceuticals, Inc

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

# combo-efficacy-safety
Supporting MATLAB code for CCR submitted manuscript: Dose Optimization for Anticancer Drug Combinations: Maximizing Therapeutic Index via Simultaneous Clinical Exposure-Toxicity/Preclinical Efficacy modeling

do_everything.m
Does almost everything. This is a script file that generates synthetic mouse efficacy and clinical toxicity data and then performs appropriate analyses with these data corresponding to some of the manuscript figures. This file calls fit_loewe_combo_tox_data2.m and unilogifit.m functions, also included. 

TGIvsGRI.m
Generates the supplemental figure comparing the TGI and GRI xenograft response metrics

fig_toxmodel.m
Generates supplementary figure that shows which alpha level corresponds to non-overlapping tox
