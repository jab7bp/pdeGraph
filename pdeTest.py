from ROOT import *
from math import sqrt
from math import log
from math import pi


# light source
# assume Poisson distribution of photons based on mean mu
# and additional Gaussian variation of mu on each pulse
class LightSource():
    def __init__(self,muGamma=20, pulseJitter=0.10):
        self.muGamma=muGamma   # average photons per pulse
        self.pulseJitter=pulseJitter # fluctuation in pulse intensity
        self.trand=TRandom(0)
        self.count=0
        self.pulsehist=TH1F()
        self.pulsehist.SetTitle("Photons per pulse")
        self.pulsehist.SetName("nPhotons")

    def Print(self):
        print "Light source parameters"
        print "<photons/pulse>",self.muGamma
        print "SD of pulse intensity:",self.pulseJitter
        
    def Pulse(self):         # pulse intensity in # of photons
        if (self.count==0):
            xmin=0
            xmax=int(self.muGamma+sqrt(self.muGamma)*4)
            self.pulsehist.SetBins(xmax,xmin,xmax)
        nPhot=self.muGamma*(1+self.trand.Gaus()*self.pulseJitter)
        nPhot=self.trand.Poisson(nPhot)
        self.pulsehist.Fill(nPhot)
        return nPhot
        
# very simple SIPM model
class SIPM():
    def __init__(self,pde=0.2,noise=0.1,enf=0.1,dcr=2,M=1,pulseWid=25):
        # device properties
        self.pde=pde
        self.noise=noise        # electronic noise relative to 1PE peak
        self.enf=enf            # excess noise factor (eg uniformity)
        self.dcr=dcr            # in MHz
        self.M=M                # Gain = counts/1pe
        self.pulseWid=pulseWid  # in ns
        self.trand=TRandom(0)

    def Print(self):
        print "---"
        print "SIPM input parameters"
        print "pde:",self.pde
        print "gain",self.M
        print "noise, enf",self.noise,self.enf
        print "---"

    def NoiseSample(self):
        return self.trand.Gaus()*self.noise*self.M

    # calculate signal contribution from DCR
    def DarkSample(self,gate):
        # dark pulses expected in integration window
        extendedWindow=self.pulseWid*1e-9*(gate+1)
        reducedWindow=self.pulseWid*1e-9*(gate-1)
        pulseAny=self.dcr*1e6 * extendedWindow
        pulseFull=self.dcr*1e6 * reducedWindow
        r=self.trand.Uniform()
        darkSig=0
        if r<pulseFull: darkSig=1
        elif r<pulseAny: darkSig=self.trand.Uniform()
        darkSig=darkSig*self.M*(1+self.trand.Gaus()*self.enf)
        return darkSig

    # simulate charge collection from pulses in gate window, including effect of random overlaps fronm DCR 
    # can also add afterpulsing, but how to model?  Probably needs recharge time + falling expo probabiity
    def SimPhD(self,pulser,gate,ntrials=20000):  # get pulse height distribution
        if gate<1: print "Warning gate < pulse width"
        muGamma=pulser.muGamma
        xmin=-3*self.noise
        xmax=self.pde*pulser.muGamma+4*sqrt(self.pde*pulser.muGamma)*self.M
        hPhD = TH1F("hPhD","Pulse height dist",300,xmin,xmax)
        for nt in range(ntrials):
            signal=self.NoiseSample()+self.DarkSample(gate) # noise+dark counts
            nPhot=pulser.Pulse()      # ignore gate assume in time, wid>=1.0
            for np in range(nPhot):
                if self.trand.Uniform()<self.pde:
                    signal=signal+self.M*(1+self.trand.Gaus()*self.enf)
            hPhD.Fill(signal)
        return hPhD

    # gate describes integration window relative to pulse width
    # eg gate = 1, means gate = pulse width
    # very simple dark count model
    def SimDarkPhD(self, gate, ntrials=10000):
        if gate<1: print "Warning gate < pulse width"
        xmin=-3*self.noise
        xmax=self.M*4
        hDark = TH1F("hDark","Dark pulse height dist",300,xmin,xmax)
        for nt in range(ntrials):
            signal=self.NoiseSample() + self.DarkSample(gate)
            hDark.Fill(signal)
        return hDark

# not used - TBD
class PeakFitter():
    def __init__(self):
        code="""
        // "simple" function to fit multiple peaks in pulse height distribution
        // par[0] : # of peaks to fit 0=noise only, 1=noise+1pe, ....
        // par[1] : noise peak normalization
        // par[2] : noise peak mean
        // par[3] : noise peak width
        // par[4] : enf
        // par[5] : gain
        // par[6] : np1 normalization
        // par[7] : np2 normalization
        // ...
        Double_t fcn(Double_t *xp, Double_t *par){
        double x=xp[0];
        int npeFit=par[0];
        double noise = par[1]*TMath::Gaus(x,par[2],par[3]);
        double val=noise;
        double enf=par[4];
        double M=par[5];
        for (int npe=1; npe<=npeFit; npe++){
          double mu=par[2]+M*npe;
          double sig=TMath::Sqrt(par[3]*par[3]+npe*enf*enf);
        val+=TMath::Gaus(x,mu,sig);
        }
        return val;
        };
        """
        gInterpreter.ProcessLine(code)

    
class PhDAnalyzier():
    def __init__(self,hPhD,hDark=None):
        self.hPhD=hPhD
        self.ts=TSpectrum()
        self.hPhD0=hDark    # dark count pulse height distributions
        self.xyPeaks=[]
        
    # use this method for the dark pulse height distribution     
    def Fit0Peak(self):
        if not self.hPhD0: return
        maxbin=self.hPhD0.GetMaximumBin()
        max=self.hPhD0.GetBinContent(maxbin)
        mu=self.hPhD0.GetBinCenter(maxbin)
        for ibin in range(maxbin-1,0,-1):
            y=self.hPhD0.GetBinContent(ibin)
            sig=mu-self.hPhD0.GetBinCenter(ibin)
            if y<max/2: break
        xmin=self.hPhD0.GetBinCenter(1)
        xmax=self.hPhD0.GetBinCenter(maxbin)+sig*1.0 # 1.0 is a hack!
        self.hPhD0.Fit("gaus","","",xmin,xmax)


    def FindPeaks(self):
        self.hfft=self.hPhD.FFT(0,"RE")
        self.hfft.SetBinContent(1,0)  # suppress DC component
        self.hfft.SetBinContent(2,0)
        self.hfft.SetBinContent(self.hfft.GetNbinsX(),0)
        self.hfft.SetBinContent(self.hfft.GetNbinsX()-1,0)
        max=self.hfft.GetMaximumBin()
        if max>self.hfft.GetNbinsX()/2: max = self.hfft.GetNbinsX()-max
        #print max
        if self.hPhD0:
            fcn=self.hPhD0.GetFunction("gaus")
            self.peakWid=fcn.GetParameter(2)/self.hPhD0.GetBinWidth(1)
            # for overlapping peaks, setting the search width smaller
            # appears to help finding the peaks
            self.peakWid=self.peakWid*0.75 
        else:  # if dark spectrum is not available
            print "Warning: dark spectrum not found"
            self.peakWid=(self.hPhD.GetNbinsX()/max)/4  # est. peak distance / 4 in Nbins
            self.peakWid=int(max/2.0)  # est. peak distance / 4 in Nbins       
        print "peakwid",self.peakWid
        #self.hfft.Draw("hist")
        #raw_input("Press Enter to continue...")
        self.npeaks=self.ts.Search(self.hPhD,self.peakWid)

        print "TSpectrum found",self.npeaks,"peaks"
        xvals=self.ts.GetPositionX()
        yvals=self.ts.GetPositionY()
        for i in range(self.npeaks):  # store peaks as a list of 2 element lists
            self.xyPeaks.append([xvals[i],yvals[i]])
        self.xyPeaks.sort()
        del self.xyPeaks[6:]   # limit analysis to first 6 peaks (including 0pe)
        # check to see if the 0PE peak was found
        if self.hPhD0:
            fcn=self.hPhD0.GetFunction("gaus")
            npe0mu=fcn.GetParameter(1)
            npe0sig=fcn.GetParameter(2)
            if abs(npe0mu-self.xyPeaks[0][0])>npe0sig:
                print "Warning: Noise peak is missing, attempting to fix..."
                xaxis=self.hPhD.GetXaxis()
                npe0bin=xaxis.FindBin(npe0mu)
                npe0y=self.hPhD.GetBinContent(npe0bin)
                print self.xyPeaks
                self.xyPeaks.insert(0,[npe0mu,npe0y])
                print self.xyPeaks
                self.npeaks=self.npeaks+1
        return self.npeaks

    # warning must call Fit0Peak and FindPeaks first
    def FitPhD(self):
        parnames=["nfit","a0","mu0","sig0","enf","gain"] # other names appended below
        gROOT.ProcessLine(".L PEFitter.C+")
        xmin=self.hPhD.GetXaxis().GetXmin()
        xmax=self.hPhD.GetXaxis().GetXmax()
        #print "fcn min max",xmin,xmax
        npefcn.SetRange(xmin,xmax)
        npefcn.SetNpx(self.hPhD.GetNbinsX())

        npePeaks=len(self.xyPeaks)-1  # remove 0pe peak from count
        fcn0=self.hPhD0.GetFunction("gaus")
        ymax=self.hPhD.GetMaximum()
        mu0=fcn0.GetParameter(1)  # noise peak mean
        sig0=fcn0.GetParameter(2)
        gain=self.xyPeaks[1][0]-self.xyPeaks[0][0] # approx gain as dist btwn peak 1 and peak 0

        npefcn.FixParameter(0,npePeaks) # of peaks to fit 0=noise only, 1=noise+1pe, ....
        npefcn.SetParameter(1,self.xyPeaks[0][1]*0) # noise peak normalization "a0"
        npefcn.SetParLimits(1,0,ymax)
        npefcn.SetParameter(2,mu0) 
        npefcn.SetParameter(3,sig0/gain) # use noise peak width as fraction of gain
        npefcn.SetParLimits(2,mu0-2*sig0,mu0+2*sig0)
        npefcn.SetParameter(4,sig0/gain) # enf -- starting guess, again as fraction of gain
        npefcn.SetParameter(5,gain) # gain approx, dist btwn peak1&0
        
        for i in range(npePeaks):
            npefcn.SetParameter(6+i,self.xyPeaks[1+i][1]) # heights of peaks 1...n
            npefcn.SetParLimits(6+i,0,ymax)
            parnames.append("a"+str(i+1))
            #print "a"+str(i+1),self.xyPeaks[1+i][1]
        for i in range(len(parnames)): npefcn.SetParName(i,parnames[i])
        for i in range(len(parnames),npefcn.GetNpar()): npefcn.FixParameter(i,0)  # unused parameters
        xend = self.xyPeaks[-1][0]+sig0*1  # end of fit range
        #c3=TCanvas()
        #npefcn.Draw()
        #raw_input("Press Enter to continue...")
        self.hPhD.Fit("npefcn","","",xmin,xend)     
        self.fcn0=TF1("fcn0","gaus",xmin,xmax)
        self.fcn0.SetParameters(npefcn.GetParameter(1),npefcn.GetParameter(2),
                                npefcn.GetParameter(3)*npefcn.GetParameter(5))
        self.fcn0.SetRange(xmin,xmax)
        self.fcn0.SetLineStyle(2)
        
    
    def CalcNpe(self):  # calculate average # of detected photons per pulse
        
        A0_w_err = []
        n0_w_err = []
        A1_w_err = []
        n1_w_err = []
        
        hPhD_test = 0
        hPhD0_Integral = 0
        
        npe_w_err = []
        # check for under/overflows
        uflow=self.hPhD.GetBinContent(0)+self.hPhD.GetBinContent(0)
        oflow=self.hPhD0.GetBinContent(self.hPhD0.GetNbinsX()+1)+self.hPhD.GetBinContent(self.hPhD.GetNbinsX()+1)
        if uflow+oflow >0: print "*** Warning, under/over flow entries in pulse height histograms"
        sigma_logA1N1 = 0
        self.Fit0Peak()
        npeaks=self.FindPeaks()
        zeropeak=(self.xyPeaks[0])[0]
        xmin=zeropeak-self.peakWid*self.hPhD.GetBinWidth(1)*1.5 # 1.5 is a hack!
        xmax=zeropeak+self.peakWid*self.hPhD.GetBinWidth(1)*1.5
        self.hPhD.Fit("gaus","","",xmin,xmax)
        fcn=self.hPhD.GetFunction("gaus")
        A=fcn.GetParameter(0)
        A0 = fcn.GetParameter(0)
        A_error = fcn.GetParError(0)
        self.A = fcn.GetParameter(0)
        self.A_err = fcn.GetParError(0)
        mu=fcn.GetParameter(1)
        mu_err = fcn.GetParError(1)
        self.noise=fcn.GetParameter(2)
        self.noise_err = fcn.GetParError(2)
        noise = fcn.GetParameter(2)
        noise0 = fcn.GetParameter(2)
        Noise_err = fcn.GetParError(2)
        
        print "Parameter A0 is: ", A0
        print "A0 Error is: ", A_error
        print "Parameter Noise0 is: ", noise0
        print "Noise0 Error is: ", Noise_err
        print "hPhD Bin Width is: ", self.hPhD.GetBinWidth(1)
        print "hPhD GetEntries is: ", self.hPhD.GetEntries()
        
        
        self.noise_err = fcn.GetParError(2)
        self.nPed=A/self.hPhD.GetBinWidth(1) * sqrt(2*pi) * self.noise
        nDarkPed=1
        nDarkTot=1
        
        if self.hPhD0:
            fcn=self.hPhD0.GetFunction("gaus")
            A=fcn.GetParameter(0)
            A1 = fcn.GetParameter(0)
            A1Err = fcn.GetParError(0)
            
            self.A1 = fcn.GetParameter(0)
            self.A1Err = fcn.GetParError(0)
            mu=fcn.GetParameter(1)
            noise=fcn.GetParameter(2)
            noise1 = fcn.GetParameter(2)
            noise1Err = fcn.GetParError(2)
            self.noise1 = fcn.GetParameter(2)
            self.noise1Err = fcn.GetParError(2)
            
            nDarkPed=A/self.hPhD0.GetBinWidth(1) * sqrt(2*pi) * noise
            nDarkTot=self.hPhD0.Integral(1,self.hPhD0.GetNbinsX()) # histogram is weighted, us integral, not # entries
            self.nDarkPed=A/self.hPhD0.GetBinWidth(1) * sqrt(2*pi) * noise
            self.nDarkTot=self.hPhD0.Integral(1,self.hPhD0.GetNbinsX())
        
            A1N1 = (self.A1*self.noise1*sqrt(2*pi))/(self.hPhD0.GetBinWidth(1)*self.hPhD0.Integral(1,self.hPhD0.GetNbinsX()))
            sigma_A1N1 = A1N1*sqrt(((self.A1Err/self.A1)**2) + ((self.noise1Err/self.noise1)**2))
            
            logA1N1 = log(A1N1)
            sigma_logA1N1 = sigma_A1N1/abs(A1N1)
            
            hPhD0_Integral = self.nDarkTot
            
            global hPhD0_Integral
            
            hPhD_test = 1
        
            print "Parameter A1 is: ", A1
            print "A1 Error is: ", A1Err
            print "Parameter Noise 1: ", noise1
            print "Noise 1 Error: ", noise1Err
            print "hPhD0 Bin Width is:", self.hPhD0.GetBinWidth(1)
            print "hPhD0 calc Integral is: ", hPhD0_Integral
        
        hPhD_BinWidth = self.hPhD.GetBinWidth(1)
        hPhD_GetEntries = self.hPhD.GetEntries()
        hPhD0_BinWidth = self.hPhD0.GetBinWidth(1)
        
        
        self.npe = -log(self.nPed/self.hPhD.GetEntries()) + log(nDarkPed/nDarkTot)
        
        AN = (self.A*self.noise*sqrt(2*pi))/(self.hPhD.GetBinWidth(1)*self.hPhD.GetEntries())
        
        sigma_AN = AN*sqrt(((self.A_err/self.A)**2) + ((self.noise_err/self.noise)**2))
        
        logAN = log(AN)
        sigma_logAN = sigma_AN/abs(AN)
        
        NPE_sigma = sqrt(((sigma_logAN)**2) + ((sigma_logA1N1)**2))
        NPE = self.npe
        
        npe_w_err.append(NPE)
        npe_w_err.append(NPE_sigma)
        
        ####Send hPhD0_Integral out as npe_w_err[2] in order to send out updated value
        npe_w_err.append(hPhD0_Integral)
        
        A0_w_err.append(A0)
        A0_w_err.append(A_error)
        
        n0_w_err.append(noise0)
        n0_w_err.append(Noise_err)
        if self.hPhD0:
            A1_w_err.append(A1)
            A1_w_err.append(A1Err)
        
            n1_w_err.append(noise1)
            n1_w_err.append(noise1Err)
        else:
            A1_w_err.append(0)
            A1_w_err.append(0)
            n1_w_err.append(0)
            n1_w_err.append(0)
            hPhD0_BinWidth = 0
        
        global A0_w_err
        global n0_w_err
        global A1_w_err
        global n1_w_err
        global hPhD_BinWidth
        global hPhD_GetEntries
        global hPhD0_BinWidth
        global hPhD_test
        
        global A
        global noise
        global NPE_sigma
            #return self.npe
        return npe_w_err

    # A fairly trivial method.  All the work is done above
    # Get parameter and their errors
    
    
    def GetPDE(self,muGamma):
        return self.npe/muGamma

    def GetNoise(self):
        noise_w_err = []
        noise_w_err.append(npefcn.GetParameter("sig0"))
        noise_w_err.append(npefcn.GetParError(3))
        return noise_w_err
                           
    def GetNoiseErr(self):
        return npefcn.GetParError(3)
    
    def GetENF(self):
        enf_w_err = []
        enf_w_err.append(npefcn.GetParameter("enf"))
        enf_w_err.append(npefcn.GetParError(4))
        return enf_w_err
    def GetENFErr(self):
        return npefcn.GetParError(4)
    
    def GetGain(self):
        gain_w_err =[]
        gain_w_err.append(npefcn.GetParameter("gain"))
        gain_w_err.append(npefcn.GetParError(5))
        #return npefcn.GetParameter("gain")
        return gain_w_err
    
    def GetGainErr(self):
        return npefcn.GetParError(5)
    def GetA(self):
        return A
    def GetA0(self):
        return A0_w_err
    def GetN0(self):
        return n0_w_err
    def GetA1(self):
        return A1_w_err
    def GetN1(self):
        return n1_w_err
    def GetNPEsigma(self):
        return NPE_sigma
    def GethPhD_BinWidth(self):
        return hPhD_BinWidth
    def GethPhD_GetEntries(self):
        return hPhD_GetEntries
    def GethPhD0_BinWidth(self):
        return hPhD0_BinWidth
    def GethPhD0_Integral(self):
        return hPhD0_Integral
    def GethPhD_test(self):
        return hPhD_test


    

if __name__ == "__main__":
    tf_out=TFile("pdeTest.root","recreate")
    pulser=LightSource()
    s=SIPM()
    s.dcr=2
    gate=1.5  # gate width relative to pulse width
    hLight=s.SimPhD(pulser,gate)
    hDark=s.SimDarkPhD(gate)
    ana=PhDAnalyzier(hLight,hDark)
    npe=ana.CalcNpe() # depends only on the ratios of events in the 0 peak
    ana.FitPhD()      # do a nice fit to the peaks


    screenY=TGClient.Instance().GetDisplayHeight()
    c1=TCanvas("results","results",int(screenY*.75),int(screenY*.75))
    c1.Divide(2,2)
    c1.cd(1)
    ana.hPhD.Draw()
    c1.cd(2)
    pulser.pulsehist.Draw()
    c1.cd(3)
    if ana.hPhD0: ana.hPhD0.Draw()
    else: hDark.Draw()
    c1.cd(4)
    #ana.FitPhD()
    npefcn.Draw()


    s.Print()

    print "Calculated PDE =",'{0:.1f}%'.format(ana.GetPDE(pulser.muGamma)*100)
    print "Calculated noise =",'{0:.2f}'.format(ana.GetNoise())
    print "Calculated ENF =",'{0:.2f}'.format(ana.GetENF())
    print "Calculated Gain =",'{0:.2f}'.format(ana.GetGain())
    

    tf_out.Write()
    tf_out.Close()
    raw_input("Press Enter to continue...")
