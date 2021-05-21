import numpy as np
import matplotlib.pyplot as plt

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
        self.sign = self.value/np.abs(self.value) if self.value != 0 else 1
        
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
    
    def __int__(self):
        return int(self.value)
    
    def __len__(self):
        return 1
    
    def __float__(self):
        return float(self.value)
    
    def __neg__(self):
        return AsymmetricUncertainty(-self.value,self.minus,self.plus)
        
    def __add__(self,other):
        if isinstance(other,type(self)):
            pass
        else:
            other = AsymmetricUncertainty(other,0,0)
        result = self.value + other.value
        pos = np.sqrt(self.plus**2 + other.plus**2)
        neg = np.sqrt(self.minus**2 + other.minus**2)
        #print("added",self,"+",other,"=",AsymmetricUncertainty(result,pos,neg))
        return AsymmetricUncertainty(result,pos,neg)
    
    def __radd__(self,other):
        if isinstance(other,type(self)):
            pass
        else:
            other = AsymmetricUncertainty(other,0,0)
        result = self.value + other.value
        pos = np.sqrt(self.plus**2 + other.plus**2)
        neg = np.sqrt(self.minus**2 + other.minus**2)
        #print("added",other,"+",self,"=",AsymmetricUncertainty(result,pos,neg))
        return AsymmetricUncertainty(result,pos,neg)
    
    def __sub__(self,other):
        if isinstance(other,type(self)):
            pass
        else:
            other = AsymmetricUncertainty(other,0,0)
        result = self.value - other.value
        pos = np.sqrt(self.plus**2 + other.minus**2)
        neg = np.sqrt(self.minus**2 + other.plus**2)
        #print("subtracted",other,"from",self,"=",AsymmetricUncertainty(result,pos,neg))
        return AsymmetricUncertainty(result,pos,neg)
    
    def __rsub__(self,other):
        if isinstance(other,type(self)):
            pass
        else:
            other = AsymmetricUncertainty(other,0,0)
        result = other.value - self.value
        pos = np.sqrt(self.minus**2 + other.plus**2)
        neg = np.sqrt(self.plus**2 + other.minus**2)
        #print("subtracted",self,"from",other,"=",AsymmetricUncertainty(result,pos,neg))
        return AsymmetricUncertainty(result,pos,neg)
    
    def __mul__(self,other):
        if isinstance(other,type(self)):
            pass
        else:
            other = AsymmetricUncertainty(other,0,0)
        
        result = self.value * other.value
        pos = np.sqrt((self.plus/self.value)**2 + (other.plus/other.value)**2) * np.abs(result)
        neg = np.sqrt((self.minus/self.value)**2 + (other.minus/other.value)**2) * np.abs(result)
        #print("multiplied",self,"by",other,"=",AsymmetricUncertainty(result,pos,neg))
        return AsymmetricUncertainty(result,pos,neg)
    
    def __rmul__(self,other):
        if isinstance(other,type(self)):
            pass
        else:
            other = AsymmetricUncertainty(other,0,0)
        
        result = self.value * other.value
        pos = np.sqrt((self.plus/self.value)**2 + (other.plus/other.value)**2) * np.abs(result)
        neg = np.sqrt((self.minus/self.value)**2 + (other.minus/other.value)**2) * np.abs(result)
        #print("multiplied",other,"by",self,"=",AsymmetricUncertainty(result,pos,neg))        
        return AsymmetricUncertainty(result,pos,neg)
    
    def __truediv__(self,other): # self divided by something
        if isinstance(other,type(self)):
            pass
        else:
            other = AsymmetricUncertainty(other,0,0)
        result = self.value / other.value
        pos = np.sqrt((self.plus/self.value)**2 + (other.minus/other.value)**2) * np.abs(result)
        neg = np.sqrt((self.minus/self.value)**2 + (other.plus/other.value)**2) * np.abs(result)
        #print("divided",self,"by",other,"=",AsymmetricUncertainty(result,pos,neg))        
        return AsymmetricUncertainty(result,pos,neg)
    
    def __rtruediv__(self,other): # something divided by self
        if isinstance(other,type(self)):
            pass
        else:
            other = AsymmetricUncertainty(other,0,0)
        result = other.value / self.value
        pos = np.sqrt((other.plus/other.value)**2 + (self.minus/self.value)**2) * np.abs(result)
        neg = np.sqrt((other.minus/other.value)**2 + (self.plus/self.value)**2) * np.abs(result)
        #print("divided",other,"by",self,"=",AsymmetricUncertainty(result,pos,neg))        
        return AsymmetricUncertainty(result,pos,neg)
    
    def __pow__(self,other): # self to the something power
        if isinstance(other,type(self)):
            pass
        else:
            other = AsymmetricUncertainty(other,0,0)
        result = self.value**other.value
        pos = np.abs(result)*np.sqrt((self.plus*other.value/self.value)**2 + (other.plus*np.log(self.value))**2)
        neg = np.abs(result)*np.sqrt((self.minus*other.value/self.value)**2 + (other.minus*np.log(self.value))**2)
        #print("raised",self,"to",other,"=",AsymmetricUncertainty(result,pos,neg))        
        return AsymmetricUncertainty(result,pos,neg)
    
    def __rpow__(self,other): # something to the self power
        if isinstance(other,type(self)):
            pass
        else:
            other = AsymmetricUncertainty(other,0,0)
        result = other.value**self.value
        pos = np.abs(result)*np.sqrt((other.plus*self.value/other.value)**2 + (self.plus*np.log(other.value))**2)
        neg = np.abs(result)*np.sqrt((other.minus*self.value/other.value)**2 + (self.minus*np.log(other.value))**2)
        #print("raised",other,"to",self,"=",AsymmetricUncertainty(result,pos,neg))                
        return AsymmetricUncertainty(result,pos,neg)
    
    def log10(self):
        result = np.log10(self.value)
        pos = self.plus/(self.value*np.log(10))
        neg = self.minus/(self.value*np.log(10))
        #print("logged",self,"=",AsymmetricUncertainty(result,pos,neg))
        return AsymmetricUncertainty(result,pos,neg)
    
    def __gt__(self,other):
        if isinstance(other,type(self)):
            pass
        else:
            other = AsymmetricUncertainty(other,0,0)
        return self.value > other.value
    
    def __lt__(self,other):
        if isinstance(other,type(self)):
            pass
        else:
            other = AsymmetricUncertainty(other,0,0)
        return self.value < other.value
    
    def __lshift__(self,other): # overloaded <<; definitively less than
        if isinstance(other,type(self)):
            pass
        else:
            other = AsymmetricUncertainty(other,0,0)
        return self.maximum < other.minimum
    
    def __rshift__(self,other): # overloaded >>; definitively greater than
        if isinstance(other,type(self)):
            pass
        else:
            other = AsymmetricUncertainty(other,0,0)
        return self.minimum > other.maximum
    
    def __le__(self,other):
        if isinstance(other,type(self)):
            pass
        else:
            other = AsymmetricUncertainty(other,0,0)
        return self.value <= other.value
    
    def __ge__(self,other):
        if isinstance(other,type(self)):
            pass
        else:
            other = AsymmetricUncertainty(other,0,0)
        return self.value >= other.value

class UncertaintyArray():
    
    def __init__(self,array=[]):
        self.as_list = list(array)
        self.as_numpy = np.array(array)
        self.flattened = self.as_numpy.flatten()
        self.shape = np.shape(self.as_numpy)
        self.ndim = self.as_numpy.ndim
        
        for i in range(len(self.flattened)):
            if isinstance(self.flattened[i], AsymmetricUncertainty):
                pass
            else:
                self.flattened[i] = AsymmetricUncertainty(self.flattened[i],0,0)
                
        self.values = [elem.value for elem in self.flattened]
                
        def __len__(self):
            return len(self.as_list)
        
        def __getitem__(self, key):
            return self.as_numpy[key]
        
        def __str__(self):
            return "foo"#str([elem for elem in self.as_list])
