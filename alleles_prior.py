"""
The MIT License

Copyright (c) 2015
The University of Texas MD Anderson Cancer Center
Hamim Zafar and Ken Chen (kchen3@mdanderson.org)

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

"""


class allele_prior:
    # refers to amplification errors (?). If so I have modified from the paper as:
    """ prior_matrix[('RS', 'Q')] represents the probability of Q being the ‘intermediate allele’ given 
    the genotype g = RS. Q is a variable that takes value from {A,T,G,C}. The term ‘intermediate allele’
    refers to the allele which is called after amplification. In the absence of any amplification errors,
    Q should be either R or S. Due to the errors introduced during preparation of the sample, Q can
    differ from both R and S. In the context of single cell sequencing data, Q accounts for the FP errors
    introduced during the amplification process"""

    def __init__(self, p):
        self.prior_matrix = {}
        self.prior_matrix[('AA', 'T')] = p
        self.prior_matrix[('AA', 'G')] = p
        self.prior_matrix[('AA', 'C')] = p
        self.prior_matrix[('AA', 'A')] = 1 - 3 * p
        self.prior_matrix[('TT', 'T')] = 1 - 3 * p
        self.prior_matrix[('TT', 'G')] = p
        self.prior_matrix[('TT', 'C')] = p
        self.prior_matrix[('TT', 'A')] = p
        self.prior_matrix[('GG', 'T')] = p
        self.prior_matrix[('GG', 'A')] = p
        self.prior_matrix[('GG', 'C')] = p
        self.prior_matrix[('GG', 'G')] = 1 - 3 * p
        self.prior_matrix[('CC', 'C')] = 1 - 3 * p
        self.prior_matrix[('CC', 'G')] = p
        self.prior_matrix[('CC', 'T')] = p
        self.prior_matrix[('CC', 'A')] = p
        self.prior_matrix[('AC', 'T')] = p
        self.prior_matrix[('AC', 'G')] = p
        self.prior_matrix[('AC', 'C')] = (1 - 2 * p) / 2
        self.prior_matrix[('AC', 'A')] = (1 - 2 * p) / 2
        self.prior_matrix[('AG', 'T')] = p
        self.prior_matrix[('AG', 'C')] = p
        self.prior_matrix[('AG', 'G')] = (1 - 2 * p) / 2
        self.prior_matrix[('AG', 'A')] = (1 - 2 * p) / 2
        self.prior_matrix[('AT', 'C')] = p
        self.prior_matrix[('AT', 'G')] = p
        self.prior_matrix[('AT', 'T')] = (1 - 2 * p) / 2
        self.prior_matrix[('AT', 'A')] = (1 - 2 * p) / 2
        self.prior_matrix[('CG', 'T')] = p
        self.prior_matrix[('CG', 'A')] = p
        self.prior_matrix[('CG', 'C')] = (1 - 2 * p) / 2
        self.prior_matrix[('CG', 'G')] = (1 - 2 * p) / 2
        self.prior_matrix[('CT', 'A')] = p
        self.prior_matrix[('CT', 'G')] = p
        self.prior_matrix[('CT', 'T')] = (1 - 2 * p) / 2
        self.prior_matrix[('CT', 'C')] = (1 - 2 * p) / 2
        self.prior_matrix[('GT', 'C')] = p
        self.prior_matrix[('GT', 'A')] = p
        self.prior_matrix[('GT', 'T')] = (1 - 2 * p) / 2
        self.prior_matrix[('GT', 'G')] = (1 - 2 * p) / 2

    def getValue(self, key):
        return self.prior_matrix[key]

    def PrintPrior(self):
        print self.prior_matrix
