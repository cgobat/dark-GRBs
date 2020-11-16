import numpy as np
import matplotlib.pyplot as plt


class custom_iter: # custom iterator class that allows for retrieval of current element w/out advancing
    def __init__(self, iterable):
        self.iterator = iter(iterable)
        self.current = None
    def __next__(self):
        try:
            self.current = next(self.iterator)
        except StopIteration:
            self.current = None
        finally:
            return self.current


class AsymmetricUncertainty:
    """
    Class for handling propagation of asymmetric uncertainties assuming a pseudo-Gaussian probability distribution
    where the plus and minus errors in each direction of the nominal value are like modified 1-sigma standard devations.
    See 2008ApJ...672..433S for more.
    
    Author: Caden Gobat, George Washington University

    Parameters
    ----------
    nominal : numeric
        the nominal value of the represented quantity
    pos_err : numeric
        the plus error on the value, as in value = nominal (+pos_err, -neg_err)
    neg_err : numeric
        the minus error on the value, as in value = nominal (+pos_err, -neg_err)

    Methods
    -------
    pdf()
        plots the probability distribution function
    
    cdf()
        plots the cumulative distribution function

    Raises
    ------
    ValueError
        if you don't pass numeric values for the arguments
    TypeError
        if you don't pass numeric values for the arguments
    """
    
    def __init__(self, nominal, pos_err, neg_err):
        self.value = float(nominal)
        self.plus = float(pos_err)
        self.minus = float(neg_err)
        self.maximum = float(nominal)+float(pos_err)
        self.minimum = float(nominal)-float(neg_err)
        
    def __str__(self):
        if np.isclose(self.plus, self.minus):
            return "{} ± {}".format(self.value,self.plus)
        else:
            return "{} (+{}, -{})".format(self.value,self.plus,self.minus)
        
    def _repr_latex_(self):
        if np.isclose(self.plus, self.minus):
            return "$%f \pm %f$" %(self.value,self.plus)
        else:
            return "$%f_{-%f}^{+%f}$" %(self.value,self.minus,self.plus)
        
    def pdf(self,num_sigma=5,discretization=100,**kwargs):
        neg_x = np.linspace(self.value-(num_sigma*self.minus),self.value,discretization)
        pos_x = np.linspace(self.value,self.value+(num_sigma*self.minus),discretization)
        p_neg = np.sqrt(2)/np.sqrt(np.pi)/(self.plus+self.minus) * np.exp(-1*(neg_x-self.value)**2 / (2*self.minus**2))
        p_pos = np.sqrt(2)/np.sqrt(np.pi)/(self.plus+self.minus) * np.exp(-1*(pos_x-self.value)**2 / (2*self.plus**2))
        plt.plot(list(neg_x)+list(pos_x),list(p_neg)+list(p_pos),**kwargs)
        plt.show()
        
    def cdf(self,num_sigma=5,discretization=100,**kwargs):
        neg_x = np.linspace(self.value-(num_sigma*self.minus),self.value,discretization)
        pos_x = np.linspace(self.value,self.value+(num_sigma*self.minus),discretization)
        p_neg = np.sqrt(2)/np.sqrt(np.pi)/(self.plus+self.minus) * np.exp(-1*(neg_x-self.value)**2 / (2*self.minus**2))
        p_pos = np.sqrt(2)/np.sqrt(np.pi)/(self.plus+self.minus) * np.exp(-1*(pos_x-self.value)**2 / (2*self.plus**2))
        x = np.array(list(neg_x)+list(pos_x))
        pdf = np.array(list(p_neg)+list(p_pos))
        cdf = np.cumsum(pdf)/np.sum(pdf)
        plt.plot(x,cdf,**kwargs)
        plt.show()
        
    def __add__(self,other):
        if isinstance(other,AsymmetricUncertainty):
            pass
        else:
            other = AsymmetricUncertainty(other,0,0)
        result = self.value + other.value
        pos = np.sqrt(self.plus**2 + other.plus**2)
        neg = np.sqrt(self.minus**2 + other.minus**2)
        return AsymmetricUncertainty(result,pos,neg)
    
    def __sub__(self,other):
        if isinstance(other,AsymmetricUncertainty):
            pass
        else:
            other = AsymmetricUncertainty(other,0,0)
        result = self.value - other.value
        pos = np.sqrt(self.plus**2 + other.minus**2)
        neg = np.sqrt(self.minus**2 + other.plus**2)
        return AsymmetricUncertainty(result,pos,neg)
    
    def __mul__(self,other):
        if isinstance(other,AsymmetricUncertainty):
            pass
        else:
            other = AsymmetricUncertainty(other,0,0)
        maximum = (self.value+self.plus) * (other.value+other.plus)
        minimum = (self.value-self.minus) * (other.value-other.minus)
        result = self.value * other.value
        pos = maximum-result
        neg = result-minimum
        return AsymmetricUncertainty(result,pos,neg)
    
    def __truediv__(self,other):
        if isinstance(other,AsymmetricUncertainty):
            pass
        else:
            other = AsymmetricUncertainty(other,0,0)
        maximum = (self.value+self.plus) / (other.value-other.minus)
        minimum = (self.value-self.minus) / (other.value+other.plus)
        result = self.value / other.value
        pos = maximum-result
        neg = result-minimum
        return AsymmetricUncertainty(result,pos,neg)
    
    def __pow__(self,other):
        if isinstance(other,AsymmetricUncertainty):
            pass
        else:
            other = AsymmetricUncertainty(other,0,0)
        return self.value**other.value
    
    def __rpow__(self,other):
        if isinstance(other,AsymmetricUncertainty):
            pass
        else:
            other = AsymmetricUncertainty(other,0,0)
        return other.value**self.value