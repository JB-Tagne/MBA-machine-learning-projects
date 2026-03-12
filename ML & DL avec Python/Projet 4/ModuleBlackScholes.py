# -*- coding: utf-8 -*-
"""
Created on Fri Dec  6 19:04:21 2019

@author: SUN
"""


####################################################
# Black Scholes pour Devise EURUSD (Garman Khlhagen)
####################################################

from math import sqrt #Importation de la fonction sqrt du module math
from math import log #Importation de la fonction log du module math
from math import exp #Importation de la fonction exp du module math

import scipy.stats as stats

#----------------------------------------
# Formules Black Scholes Prime Call & Put
#----------------------------------------

# Calcul de d1 et d2
#-------------------

def d_un(S,K,sigma,T,r,q):
    d1=(log(S/K)+(r-q+(sigma**2)/2)*T) / (sigma*sqrt(T))
    return(d1)

def d_deux(S,K,sigma,T,r,q):
    d2 =d_un(S,K,sigma,T,r,q) - (sigma*sqrt(T))
    return(d2)

# Fonction de répartition N() = stats.norm.cdf(x, 0, 1)

def N(x):
    return stats.norm.cdf(x, 0, 1)

# Option Call
#------------

def bs_call(S,K,sigma,T,r,q):
    d1=d_un(S,K,sigma,T,r,q)
    d2=d_deux(S,K,sigma,T,r,q)
    
    op_call = (S*exp(-q*T)*N(d1)) - (K*exp(-r*T)*N(d2))
    
    return(op_call)

# Option put
#-----------

def bs_put(S,K,sigma,T,r,q):
    d1=d_un(S,K,sigma,T,r,q)
    d2=d_deux(S,K,sigma,T,r,q)
    
    op_put=(-S*exp(-q*T)*N(-d1))+(K*exp(-r*T)*N(-d2))
    
    return(op_put)

# Test
# ----

if __name__ == "__main__":
    
    S = 1.135 # Spot EURUSD
    r = 2.60/100 # Taux sans risque $
    q = -0.31/100 # Taux sans risque €
    T = 0.25 # 3 mois

    K = 1.1432
    sigma = 6.16/100

    prime_call = bs_call(S,K,sigma,T,r,q)
    print("prime_call:" "%.5f" % prime_call)

    prime_put = bs_put(S,K,sigma,T,r,q)
    print("prime_put:" "%.5f" % prime_put)
