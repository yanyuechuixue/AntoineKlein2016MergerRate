import os,sys,random
import numpy as np
from astropy.cosmology import Planck18_arXiv_v2 as cosmos

# codes form https://bitbucket.org/radboudradiolab/gwtoolbox
# Data from https://bitbucket.org/radboudradiolab/gwtoolbox/src/master/gwtoolbox/catalogues_mHz/MBHs/
# Paper is 10.1103/PhysRevD.93.024003

MBHBunits = {'EclipticLatitude':                 'Radian',\
         'EclipticLongitude':                'Radian',\
         'PolarAngleOfSpin1':                'Radian',\
         'PolarAngleOfSpin2':                'Radian',\
         'AzimuthalAngleOfSpin1':            'Radian',\
         'AzimuthalAngleOfSpin2':            'Radian',\
         'Spin1':                            'MassSquared',\
         'Spin2':                            'MassSquared',\
         'Mass1':                            'SolarMass',\
         'Mass2':                            'SolarMass',\
         'CoalescenceTime':                  'Second',\
         'PhaseAtCoalescence':               'Radian',\
         'InitialPolarAngleL':               'Radian',\
         'InitialAzimuthalAngleL':           'Radian',\
         'Approximant':                                      'ModelName',\
         'Cadence':                          'Seconds',\
         'Redshift':                         'dimensionless',\
         'Distance':                         'Gpc',
         'ObservationDuration':              'Seconds',
         'SourceType':                       'name'}

def Str(a):
    if type(a)==bytes:
        return a.decode('utf8')
    if type(a)!=str:
        return str(a)   
    return a

class ParsUnits():
    """
    This class defines the small object to manage parameters and units.
    """

    def __init__(self, pars_i=None, units_i=None, name='', value=0., unit=''):
        """
        Initalize `ParsUnits`
        @param pars_i is an optional dictionary of parameters
        @param units_i is an optional dictionary of units (same size as pars_i)
        @param name is an optional name
        @param value is an optional value
        @param unit is an optional unit
        """
        self.pars = {}
        self.units = {}
        if pars_i is not None:
            self.addDict(pars_i,units_i)
        elif name!='':
            self.addPar(name,value,unit)
    def copy(self, p):
        self.pars=p.pars.copy()
        self.units=p.units.copy()
        pass
    def __del__(self):
        # so far is empty, see what we need to add here
        pass

    def display(self,ReturnStr=False):
        """
        Display all parameters
        @param ReturnStr is true to return string instead of display
        """
        r = ""
        for i,k in enumerate(self.pars):
            r = r + "\t"+str(k)+" "
            if type(self.pars[k])==str:
                r = r + self.pars[k]
            else:
                r = r + str(self.pars[k])
            r = r + " ["+Str(self.units[k])+"]\n"
        if ReturnStr:
            return r
        else:
            print(r)

    def addPar(self,name,value,unit):
        """
        Add parameter, its value and unit
        @param name is the name of the parameter
        @param value is the value of the parameter
        @param unit is the unit of the parameter
        """
        self.pars.update({name : value})
        self.units.update({name : unit})

    def addDict(self,pars_i,units_i):
        """
        Add dictionnary and its unit
        @param pars_i is a dictionary of parameters
        @param units_i is a dictionary of units(same size as pars_i)
        """
        if len(pars_i)==len(units_i):
            self.pars.update(pars_i)
            self.units.update(units_i)
        else:
            raise Exception('addDict : parameters and units should have the same number of elements.')

    def get(self,parName):
        """
        Get parameter value
        @param parName parameter name
        @return value
        """
        ### TODO : Use proper error system
        if parName not in self.pars :
            print("WARNING: ParsUnits.get :",parName,"is not a parameter name (",self.pars,")")
            return None
        else:
            return self.pars[parName]

    def getConvert(self,parName,conversion,requiredUnit):
        """
        Get parameter value for distance parameters
        @param parName is the parameter name
        @param conversion is the conversion dictionary:
            + LC.convT for mass, time and distance (everything in time)
            + LC.convMass for mass only
            + LC.convTime for time only
            + LC.convDistance for distance only
        @param requiredUnit is the required unit [default:s]
        @return value
        """
        ### TODO : Use proper error system
        v = self.get(parName)
        if type(v)!=type(np.zeros(10)) and v == None:
            return None
        else:
            uV = self.units[parName]
            uV = uV.lower()
            requiredUnit = requiredUnit.lower()
            if requiredUnit not in conversion:
                print("WARNING: ParsUnits.getConvert : parameter unit",requiredUnit,"is not in the conversion list (",conversion,")")
                return None
            if uV not in conversion:
                print("WARNING: ParsUnits.getConvert : required unit",uV,"is not in the conversion list (",conversion,")")
                return None
            return v * ( conversion[uV] / conversion[requiredUnit] )
        

def readfromcat(cosmos=None, SourceType='MBHB', Approximant='IMRPhenomD', Cadence=10, duration=1, file_name=None, n=0):
    """
read event from catalgoue, and return a p object
@cosmos: an astropy.cosmo object 
@param file_name: filename of catalogue
@param n: number of event reading from the catalogue file    
@param duration: observation duration in unit of year, defalt is 1 year. 
return a p object, with 20 parameters in it. 
    """
#    masscut=1e4 
    Cat={}   
    data_cube=np.loadtxt(file_name, dtype=float, max_rows=n)
#    print("hello")
#    print('filename',file_name)
    zs=data_cube[:,0]
    Cat['SourceType']=SourceType       # 1
    Cat['Approximant']=Approximant     # 2
    Cat['Cadence']=Cadence             # 3
    Cat['Redshift']=zs                 # 4
    #mass1=(1.0+zs)*data_cube[:,1] # read in intrinsic masses, convert to red-shifted: may be need not # yishuxu: 28 Apr
    #mass2=(1.0+zs)*data_cube[:,2]
    mass1=data_cube[:,1]*(1.0+zs)
    mass2=data_cube[:,2]*(1.0+zs)
    #indices= np.array(mass1)>1e4 and np.array(mass2)>1e4
    #Cat['Mass1']=mass1[(mass1>masscut)*(mass2>masscut)]                 # 5
    Cat['Mass1']=mass1
    Cat['Mass2']=mass2 #[(mass1>masscut)*(mass2>masscut)]                 # 6
    spin1=data_cube[:,3]    
    Cat['Spin1']=spin1 #[(mass1>masscut)*(mass2>masscut)]                 # 7
    spin2=data_cube[:,4]
    Cat['Spin2']=spin2 #[(mass1>masscut)*(mass2>masscut)]                 # 8
    PolarAngleOfSpin1=data_cube[:,5]   
    PolarAngleOfSpin2=data_cube[:,6]  
    Cat['PolarAngleOfSpin1']=PolarAngleOfSpin1 #[(mass1>masscut)*(mass2>masscut)] # 9
    Cat['PolarAngleOfSpin2']=PolarAngleOfSpin2 #[(mass1>masscut)*(mass2>masscut)] #10
    Sphi12s=data_cube[:,7]
    CoalescenceTime=np.random.uniform(low=0,high=duration,size=len(mass1))*31556926. # want seconds
    Cat['CoalescenceTime']= CoalescenceTime*0.5 #[(mass1>masscut)*(mass2>masscut)] #11
    EclipticLatitude=0.5*np.pi - data_cube[:, 9] 
    Cat['EclipticLatitude']=EclipticLatitude #[(mass1>masscut)*(mass2>masscut)] #12
    EclipticLongitude=data_cube[:,10] 
    Cat['EclipticLongitude']=EclipticLongitude #[(mass1>masscut)*(mass2>masscut)] #13
    InitialPolarAngleL=data_cube[:,11]
    Cat['InitialPolarAngleL']=InitialPolarAngleL #[(mass1>masscut)*(mass2>masscut)] #14
    ObservationDuration=CoalescenceTime
    Cat['ObservationDuration']=ObservationDuration #[(mass1>masscut)*(mass2>masscut)] #15
    AzimuthalAngleOfSpin1, AzimuthalAngleOfSpin2, InitialAzimuthalAngleL, PhaseAtCoalescence =  np.random.uniform(low=0.0, high=2.*np.pi, size=(4,len(mass1)))
    Cat['AzimuthalAngleOfSpin1']=AzimuthalAngleOfSpin1 #[(mass1>masscut)*(mass2>masscut)] # 16
    Cat['AzimuthalAngleOfSpin2']=AzimuthalAngleOfSpin2 #[(mass1>masscut)*(mass2>masscut)] # 17
    Cat['InitialAzimuthalAngleL']=InitialAzimuthalAngleL #[(mass1>masscut)*(mass2>masscut)] #18
    Cat['PhaseAtCoalescence']=PhaseAtCoalescence #[(mass1>masscut)*(mass2>masscut)]  #19
    Cat['Distance']=[cosmos.luminosity_distance(z).value/1e3 for z in zs] #[(mass1>masscut)*(mass2>masscut)]]
    ps=[]

    for i in range(0,min(n,len(data_cube))): # [(mass1>masscut)*(mass2>masscut)]))):
       # ParsValues=MBHBunits
        ParsValues={}
        for kys in MBHBunits.keys():
            if type(Cat[kys])==np.ndarray or type(Cat[kys])==list:
                ParsValues[kys]=Cat[kys][i]
            else: 
                #print(i)
                ParsValues[kys]=Cat[kys]
        p=ParsUnits(pars_i=ParsValues, units_i=MBHBunits)
        ps.append(p)
    return ps

def EventsUniverse(cosmos=cosmos, duration=None, sourceType='MBHB', model=None):
    '''
    @cosmos : an astropy.cosmos object, typically could be Planck18_arxiv_v2;
    @duration : observation duration in unit of year, defalt is None;
    @sourceType : Could only be 'MBHB', otherwise will retun None;
    @model : string, can be 'pop3', 'Q3_delays', 'Q3_nodelays';
    return a ParsUnits object. Whose parameters could be gotten by ParsUnits.display(), ParsUnits.pars, ParsUnits.units, etc.
    '''
    # the duration needs to be rescaled, be cause the numbers in the catalogues seems 5 times larger than needed
    duration=duration/5.
    if duration<=0.1:
        raise ValueError('duration too short (>=0.1 years)!')
    elif duration>=10:
        raise ValueError('duration too long (<=10 years)!')
    if sourceType=='MBHB':
        path_sourcetype='/Users/zzc/temp/gwtoolbox/gwtoolbox/catalogues_mHz/MBHs/'
    else: return None
    path_catalogue=path_sourcetype+model+'/'
    all_names=os.listdir(path_catalogue)
    years_int=int(duration) # duration in unit of years, integer parts
    years_fraction=duration-years_int
    name_selected=random.sample(all_names, k=years_int+1)
    ps=[]
    cadence=10
    for file_ in name_selected[1:]:
        ps=ps+readfromcat(cosmos=cosmos, SourceType=sourceType, Approximant='IMRPhenomD', Cadence=cadence, duration=duration*5, file_name=path_catalogue+file_, n=1000)
    file_=name_selected[0]
    data_cube=np.loadtxt(path_catalogue+file_)
    num=int(years_fraction*len(data_cube))+2
    #print(num)
    ps=ps+readfromcat(cosmos=cosmos, SourceType=sourceType, Approximant='IMRPhenomD', Cadence=cadence, duration=duration*5, file_name=path_catalogue+file_, n=num)
    return ps