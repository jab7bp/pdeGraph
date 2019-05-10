####.ROOT FILE Nomenclature for parsing:
# EX: 09-28-18-D15-Ch12-1300dac-57_0V-1_250nA
#[0]: Month; MM
#[1]: Day; DD
#[2]: Year; YY
#[3]: Device ID: D15
#[4]: Channel Number; Ch12
#[5]: Laser Dac Setting; 1300dac
#[6]: Bias Voltage; 57_0V (57.0 volts)
#[7]: Reference Diode Current; 1_250nA (1.250 nA)   --> Can also be formatted in pA

# keep ROOT TApplication from grabbing -h flag
from ROOT import PyConfig
PyConfig.IgnoreCommandLineOptions = True
from ROOT import *

from pdeTest import *
import sys,os
import csv
import argparse
import numpy
import math


def roundup(x, n):
    return int(math.ceil(x/(1.0*(10**n)))*(1.0*(10**n)))

def roundDown(x,n):
    return math.floor(x*(1.0*(10**n)))/(1.0*(10**n))

def roundUpDec(x,n):
    return math.ceil(x*(1.0*(10**n)))/(1.0*(10**n))


def roundTo500(x):
    return math.ceil(x/500.0)*500.0

if len(sys.argv)<2:
    print "No input file/template given"
    print "Usage: pdeFitParams.py filename or template"
    sys.exit()

parser = argparse.ArgumentParser(description='pdeFitParams analyzer') 
parser.add_argument('files', nargs='*', help="Give file paths to specify files.")
parser.add_argument("-p", "--plotAll", default=None, help="show the plots",
                    action="store_true")
parser.add_argument('-c', '--calib', type=float, default=-1,
                        help="Ref. current to photons calibration")
parser.add_argument("-o", "--output", default=0, help = "Set file name for output of NPE fit parameters")
parser.add_argument("-u", "--uncert", default=0, help = "Set custom uncertainty value for NPE fit (added in quadrature)")


##Output file for fit parameters and values for calculating

args = parser.parse_args()
calib=args.calib
output = args.output

vbias=[]
refcurrent=[]
picorange=[]
npeval=[]
vex=[]
v_ov=[]
params=[]
freq=10000.0
ratio = 4633.39
resp = .427085
wlen = 660E-9
h = 6.62607E-34
c = 2.99E8
phot_e = h*c/wlen
nphot_ref=[]
nphot_test=[]
pde=[]

gain_arr=[]
gain_err_arr = []
noise_arr=[]
noise_err_arr = []
enf_arr=[]
enf_err_arr = []


########Uncertainty Values
wlen_unc = 10E-9
resp_unc = 0.01*resp
phot_e_unc = (wlen_unc/wlen)*phot_e

######Current Uncertainty must be defined below, circa line 320.
######Not sure why it won't take the value here, needs it in the if statement.
#current_unc = 0.03
freq_unc = 500.0
ratio_unc = 62.1828
fre = freq*resp*phot_e

## Vbias uncertainty: +/- 0.015% + 2.4mV
#fre_unc = fre*sqrt((((.01)**2) + ((phot_e_unc/phot_e)**2) + ((freq_unc/freq)**2)))

fre_unc = (fre)*sqrt(((freq_unc/freq)**2) + ((resp_unc/resp)**2) + ((phot_e_unc/phot_e)**2))

##########The following uncertainty is one to be added to accomodate/characterize
##########for the uncertainty due to the fitting algorithm
##########It will be added in quadrature to the error for the NPE fit/parameters
if(args.uncert):
    custom_unc = (float(args.uncert))
else:
    custom_unc = 0

refcurr_unc = []
nphot_ref_unc = []
nphot_test_unc = []
pde_unc = []
npe_unc=[]
vbias_unc = []
vex_unc = []
gain_unc = []

a0 = []
a0_err = []
n0 = []
n0_err = []
a1 = []
a1_err = []
n1 = []
n1_err = []

hPhD_BinWidth = []
hPhD_GetEntries = []
hPhD0_BinWidth = []
hPhD0_Integral = []

if (args.output):
    print "output file name is: " + args.output
    file = open(args.output, "w+")

else:
    file = open("pde_parameters.txt", "w+")

file.writelines("Frequency:                     %d\n" % freq)
file.writelines("Frequency Uncertainty:         %d\n" % freq_unc)
file.writelines("Wavelength:                    %e\n" % wlen)
file.writelines("Wavelength Uncertainty:        %e\n" % wlen_unc)
file.writelines("Responsitivity:                %f\n" % resp)
#file.writelines("PiN Current Uncertainty:       %f\n\n" % current_unc)

print "Ref. current to photons calibration",calib
tgGain=TGraphErrors()
tgGain.SetTitle("Gain vs. Vbias;Vbias [V];Gain [arbitrary units]")

tgPDE=TGraphErrors()
tgPDE.SetTitle("PDE vs V_ex; V_ex [V]; PDE [%]")

tgVbr=TGraphErrors()

tgVex=TGraphErrors()

tgNoise=TGraphErrors()

tgENF=TGraphErrors()

for fn in args.files:
    print fn
    parsed=os.path.basename(fn).replace(".root","").split("-")
    #print "parsed: ", parsed
    dac = "{:.0f}".format((float(parsed[5].replace("dac",""))))
    vbias.append(float(parsed[6].replace("_",".").replace("V","")))
    
    if "pa" in parsed[7]:
        refcurrent.append(float(parsed[7].replace("_",".").replace("pa",""))/1000)
    elif "na" in parsed[7]:
        refcurrent.append(float(parsed[7].replace("_",".").replace("na","")))
    elif "nA" in parsed[7]:
        refcurrent.append(float(parsed[7].replace("_",".").replace("nA","")))
    chip_id = parsed[3]
    chip_ch = parsed[4]

##Calculate uncertainty for Vbias from source meter
for v in range(len(vbias)):
    vbias_unc.append((vbias[v]*(.015/100)) + (2.4E-3))

for i in range(len(args.files)):
    
    file.write('**********************************************\n\n')
    file.writelines("%s \n" % args.files[i])

    tf=TFile(args.files[i])
    hLight=tf.Get("hpulses1")
    hDark=tf.Get("hpulses0")
    hRange=tf.Get("hRange")
    picorange.append(hRange.GetBinContent(1))

    if (args.plotAll):
        hLight.Draw()
        hDark.Draw("same")

    ana=PhDAnalyzier(hLight.Clone(),hDark.Clone())

###npe[0] = NPE calc; npe[1] = NPE Uncertainty; npe[2] = hPhD0_Integral
    npe=ana.CalcNpe()
    a0_list = ana.GetA0()
    n0_list = ana.GetN0()
    a1_list = ana.GetA1()
    n1_list = ana.GetN1()
    hPhD_test = ana.GethPhD_test()

    ana.FitPhD() # do a nice fit to the peaks

    if (args.plotAll):
        screenY=TGClient.Instance().GetDisplayHeight()
        c1=TCanvas("results","results",int(screenY*.75),int(screenY*.75))
        c1.Divide(1,2)
        c1.cd(1)
        ana.hPhD.Draw()
        c1.cd(2)
        ana.hPhD0.Draw()

    print "\nMean NPE detectected",npe[0], " +/- ", sqrt((npe[1]**2) + (custom_unc**2)), "\n"

    if (args.uncert) or not(custom_unc == 0):
        print "Custom uncertainty added in quadrature: ", custom_unc

###Store fitting NPE fitting parameters A0, N0, A1, N1 and associated errors.

    a0.append(a0_list[0])
    a0_err.append(a0_list[1])
    n0.append(n0_list[0])
    n0_err.append(n0_list[1])
    a1.append(a1_list[0])
    a1_err.append(a1_list[1])
    n1.append(n1_list[0])
    n1_err.append(n1_list[1])

    hPhD_BinWidth.append(ana.GethPhD_BinWidth())
    hPhD_GetEntries.append(ana.GethPhD_GetEntries())
    hPhD0_BinWidth.append(ana.GethPhD0_BinWidth())
    hPhD0_Integral.append(ana.GethPhD0_Integral())

####Calculated NPE Values
    npe_calc = -log(a0_list[0]*sqrt(2*pi)*n0_list[0]/(hPhD_BinWidth[0]*hPhD_GetEntries[0])) + log(a1_list[0]*sqrt(2*pi)*n1_list[0]/(hPhD0_BinWidth[0]*npe[2]))

###Write the parameters out to a file
    file.writelines("\nParameter A0:          %s \n" % a0_list[0])
    file.writelines("A0 Error:              %s \n" % a0_list[1])
    file.writelines("Parameter N0:          %s \n" % n0_list[0])
    file.writelines("N0 Error:              %s \n" % n0_list[1])
    file.writelines("hPhD Bin Width:        %s \n" % hPhD_BinWidth[0])
    file.writelines("hPhD Get Entries:      %s \n\n" % hPhD_GetEntries[0])

    file.writelines("Parameter A1:          %s \n" % a1_list[0])
    file.writelines("A1 Error:              %s \n" % a1_list[1])
    file.writelines("Parameter N1:          %s \n" % n1_list[0])
    file.writelines("N1 Error:              %s \n" % n1_list[1])
    file.writelines("hPhD0 Bin Width:       %s \n" % hPhD0_BinWidth[0])
    file.writelines("hPhD0 Integral:        %s \n\n" % npe[2])

    if (args.uncert) or not(custom_unc == 0):
        file.writelines("Custom uncertainty:    %s \n\n" % custom_unc)

    if (args.plotAll): raw_input("Press Enter to continue...")

    gain_raw=ana.GetGain()
    gain = gain_raw[0]*picorange[i]/100 # scale all gains to 100 mV range
    gain_err = gain_raw[1]*picorange[i]/100 # scale all gain errors to 100 mV range

    tgGain.SetPoint(i,vbias[i],gain)

    row=[sys.argv[1],npe,ana.GetNoise(),ana.GetGain(),ana.GetENF(),tf.Get("hRange").GetBinContent(1)]
    

    npeval.append(npe[0]) #Fill NPE array
    npe_unc.append(sqrt((npe[1]**2) + (custom_unc**2))) #Fill NPE uncertainty array plus Custom Unc.

    gain_arr.append(gain) #Fill Gain array
    gain_err_arr.append(gain_err)
    noise_raw = ana.GetNoise()
    noise_arr.append(noise_raw[0])
    noise_err_arr.append(noise_raw[1])
    
    enf_raw = ana.GetENF()
    enf_arr.append(enf_raw[0])
    enf_err_arr.append(enf_raw[1])


#tcan=TCanvas("fitparams","Fit params")
tgGain.Fit("pol1", "P", "", vbias[0], vbias[len(vbias)-1])
#tgGain.Draw("ALP")

func = tgGain.GetFunction("pol1")

#Determine x-intercept -> V_breakdown
p0 = func.GetParameter(0);
p0_err = func.GetParError(0);
p1 = func.GetParameter(1);
p1_err = func.GetParError(1);


vbr =  -round(p0/p1,2)
vbr_unc = vbr*sqrt(((p0_err/p0)**2) + ((p1_err/p1)**2))


ylast = round(p1*(vbias[len(vbias)-1])+p0,2)

#Define points of extrapolated line
vbr_fit = [(vbr,0), vbias[len(vbias)-1],ylast]

#vbr_x =[vbr, vbias[(len(vbias)-1)]]
vbr_x = [vbr, ((ylast - p0)/(p1))]
vbr_y =[0, ylast]


#Set points to draw extrapolated V_br line
for i in range(len(vbr_x)):
    tgVbr.SetPoint(i,vbr_x[i],vbr_y[i])

#### Draw VBR on Gain vs V_Bias graph
"""tgVbr.Draw("AlP")

tcan.Draw()"""

#Gain vs. Voltage graph settings
tgGain.GetXaxis().SetLimits(vbr*.99,1.005*vbias[len(vbias)-1])
tgGain.GetYaxis().SetRangeUser(0,roundup(max(gain_arr) + gain_err_arr[gain_arr.index(max(gain_arr))], 3))
"""tgGain.Draw("ALP")"""
tgVbr.GetXaxis().SetLimits(vbr*.99,vbias[len(vbias)-1])
tgVbr.GetYaxis().SetRangeUser(0,roundup(gain_arr[len(gain_arr)-1], 3))
"""tgVbr.Draw("same")"""
tgVbr.SetLineColor(2)
tgVbr.SetLineWidth(2)

text = TLatex()
text.SetTextAlign(13)

#### Draw VBR on Gain vs V_Bias graph
"""text.DrawLatex(0.995*vbr,-500,Form("V_br = %g +/- %g V" % (vbr,round(vbr_unc,2))))

text.DrawClone()
tcan.Modified()"""

#Determine number of photons at test/ref per each ref. curr.

for i in range(len(args.files)):
    current_unc = 0.03
    vex.append(vbias[i] - vbr)
    vex_unc.append(sqrt(((vbias_unc[i])**2) + ((vbr_unc)**2)))
    #vex_unc.append((vbr_unc/vbr)*vex[i])
    
    nphot_ref.append((refcurrent[i]*1E-9)/(resp*freq*phot_e))
    refcurr_unc.append(current_unc*refcurrent[i])
    nphot_test.append(nphot_ref[i]/ratio)

    nphot_ref_unc.append(nphot_ref[i]*sqrt(((current_unc)**2) + ((fre_unc/fre)**2)))
    nphot_test_unc.append(nphot_test[i]*sqrt(((nphot_ref_unc[i]/nphot_ref[i])**2)+ ((ratio_unc/ratio)**2)))

print "\nBreakdown voltage: ", vbr, " +/- ", round(vbr_unc, 3), " [V]"
print "Ratio: ", ratio, " +/- ", ratio_unc, "\n"

#From number of photons determine NPE at DUT
##Then plot various parameters
for i in range(len(args.files)):
    pde.append(npeval[i]/nphot_test[i])
    #set points for plotting, PDE, Gain,ENF,  and Noise v V_ex
    tgPDE.SetPoint(i,vex[i], pde[i]*100)
    tgVex.SetPoint(i,vex[i], gain_arr[i])
    tgNoise.SetPoint(i,vex[i],noise_arr[i])
    tgENF.SetPoint(i,vex[i],enf_arr[i])

    #uncertainty for PDE
    pde_unc.append(pde[i]*sqrt(((nphot_test_unc[i]/nphot_test[i])**2) + ((npe_unc[i]/npeval[i])**2)))
    

    print "PDE and uncertainty: ", pde[i], " +/- ", pde_unc[i]

###Draw PDE vs V_EX Graph
"""tcanpde=TCanvas("pde", "PDE vs V_ex")
#tgPDE.Fit("pol3")
tgPDE.Draw("ALP")
tgPDE.GetXaxis().SetLimits(roundDown(0.95*(min(vex) - vex_unc[vex.index(min(vex))]), 1), roundUpDec(1.05*(max(vex) + vex_unc[vex.index(max(vex))]),1))
tgPDE.GetYaxis().SetRangeUser(roundDown(95.0*(min(pde) - pde_unc[pde.index(min(pde))]),3), roundUpDec(105*(max(pde) + pde_unc[pde.index(max(pde))]), 3))
#tgPDE.GetYaxis().SetRangeUser(90*pde[0], 110*pde[len(pde)-1])
tcanpde.Draw()"""

#tcanvex=TCanvas("vex", "Gain vs V_ex")
#tgVex.Draw("ALP")
#tcanvex.Draw()

tgVex.GetXaxis().SetTitle("V_ex [V]")
tgVex.GetXaxis().SetLimits(roundDown(0.99*(min(vex) - vex_unc[vex.index(min(vex))]), 1), roundUpDec(1.01*(max(vex) + vex_unc[vex.index(max(vex))]),1))

tgVex.GetYaxis().SetTitle("Gain [A. U.]")
tgVex.GetYaxis().SetRangeUser(0.9*gain_arr[0], roundTo500(gain_arr[len(gain_arr)-1]))
tgVex.SetTitle("Gain vs V_ex")

##Putting all plots onto a single split canvas


allcan = TCanvas("allcan", '{0} {1} PDE Plots ({2} dac)'.format(chip_id, chip_ch, dac),2000,1700)

#Canvas Title
allcanLabel = TPaveLabel(.1,.96,.9,.99, '{0} {1} PDE Plots ({2} dac)'.format(chip_id, chip_ch,dac))
allcanLabel.Draw()
graphPad = TPad("Graphs", "Graphs",0.06,0.05,0.95,.95)
graphPad.Draw()
graphPad.cd()
graphPad.Divide(3,2)


#plot pde vs V_ex
graphPad.cd(1)
tgPDE.Draw()

#plot gain v V_ex
graphPad.cd(2)
tgVex.Draw()

#plot Gain v V_bias
graphPad.cd(3)
tgGain.Draw()
tgVbr.Draw("same")
text = TLatex()
text.SetTextAlign(13)
text.DrawLatex(0.995*vbr,-500,Form("V_br = %g +/- %g V" % (vbr, round(vbr_unc,2))))

#Plot Noise v V_ex
graphPad.cd(4)
tgNoise.Draw("ALP")
tgNoise.GetXaxis().SetTitle("V_ex [V]")
tgNoise.GetXaxis().SetLimits(roundDown(0.95*(min(vex) - vex_unc[vex.index(min(vex))]), 1), roundUpDec(1.05*(max(vex) + vex_unc[vex.index(max(vex))]),1))
tgNoise.GetYaxis().SetTitle("Noise [A. U.]")
tgNoise.GetYaxis().SetRangeUser(0.95*(min(noise_arr) - noise_err_arr[noise_arr.index(min(noise_arr))]), round(1.05*(max(noise_arr) + noise_err_arr[noise_arr.index(max(noise_arr))]),2))
tgNoise.SetTitle("Noise vs V_ex")

#Plot ENF vs V_ex
graphPad.cd(5)
tgENF.Draw("ALP")
tgENF.GetXaxis().SetTitle("V_ex [V]")
tgENF.GetXaxis().SetLimits(roundDown(0.95*(min(vex) - vex_unc[vex.index(min(vex))]), 1), roundUpDec(1.05*(max(vex) + vex_unc[vex.index(max(vex))]),1))
tgENF.GetYaxis().SetTitle("ENF [A. U]")
tgENF.GetYaxis().SetRangeUser(0.99*(min(enf_arr) - enf_err_arr[enf_arr.index(min(enf_arr))]), 1.01*(max(enf_arr) + enf_err_arr[enf_arr.index(max(enf_arr))]))
tgENF.SetTitle("ENF vs V_ex")
#graphPad.cd(5).SetLeftMargin(.15)

for i in range(1,6):
    graphPad.cd(i).SetLeftMargin(.15)



#c2 = TCanvas('c2', 'Graph with errors', 200, 10, 700, 500)

gr = TGraphErrors()

for i in range(len(pde)):
    gr.SetPoint(i, vex[i], pde[i]*100)
    gr.SetPointError(i, vex_unc[i], pde_unc[i]*100)
    tgPDE.SetPointError(i, vex_unc[i], pde_unc[i]*100)
    tgGain.SetPointError(i, vbias_unc[i], gain_err_arr[i])
    tgVex.SetPointError(i, vex_unc[i], gain_err_arr[i])
    tgNoise.SetPointError(i, vex_unc[i], noise_err_arr[i])
    tgENF.SetPointError(i, vex_unc[i], enf_err_arr[i])

#tgGain.Fit("pol1", "P")

#gr.Draw("ALP")
#c2.Draw()


raw_input("\nPress Enter to continue...")
