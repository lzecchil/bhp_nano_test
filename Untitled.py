#!/usr/bin/env python
# coding: utf-8

# In[348]:


import os
import ROOT
from ROOT import RDataFrame
from ROOT import TChain, TSelector, TTree, TH1F, TCanvas, TFile, TEfficiency, TLegend
from ROOT import TLorentzVector
from ROOT import TVector3
from array import array
import numpy as np
ROOT.gROOT.SetBatch(ROOT.kTRUE)
file ="../bph_nano_test/bph_nano_test_mc_TracksMuons.root"
f=TFile(file)
tree=f.Get("Events")
#tree.Print()
N = tree.GetEntries()
print("total number of events=",N)


# In[349]:


h_IM_Phi=TH1F("h_IM_Phi","K+K- invariant mass",100,0.98,1.1)
h_IM_Ks=TH1F("h_IM_Phi","pi+pi- invariant mass",100,0.2,1)
h_IM_B=TH1F("h_IM_B","K+K-pi+pi- invariant mass",100,4.6,6)
h_pt_Kp=TH1F("h_pt_Kp","reco K+ pt distribution",100,0,25)
h_pt_Km=TH1F("h_pt_Km","reco K- pt distribution",100,0,25)
h_pt_pip=TH1F("h_pt_pip","reco pi+ pt distribution",100,0,25)
h_pt_pim=TH1F("h_pt_pim","reco pi- pt distribution",100,0,25)
h_Eta_Kp=TH1F("h_Eta_Kp","reco K+ Eta distribution",100,-3,3)
h_Eta_Km=TH1F("h_Eta_Km","reco K- Eta distribution",100,-3,3)
h_Eta_pip=TH1F("h_Eta_pip","reco pi+ Eta distribution",100,-3,3)
h_Eta_pim=TH1F("h_Eta_pim","reco pi- Eta distribution",100,-3,3)
h_Phi_Kp=TH1F("h_Phi_Kp","reco K+ Phi distribution",100,-4,4)
h_Phi_Km=TH1F("h_Phi_Km","reco K- Phi distribution",100,-4,4)
h_Phi_pip=TH1F("h_Phi_pip","reco pi+ Phi distribution",100,-4,4)
h_Phi_pim=TH1F("h_Phi_pim","reco pi- Phi distribution",100,-4,4)
hNumK=TH1F("hNumK","reco K vs pt",100,0,25)
hNumKp=TH1F("hNumKp","reco K+ vs pt",100,0,25)
hNumKm=TH1F("hNumKm","reco K- vs pt",100,0,25)
hNumpi=TH1F("hNumpi","reco pi vs pt",100,0,25)
hNumpip=TH1F("hNumpip","reco pi+ vs pt",100,0,25)
hNumpim=TH1F("hNumpim","reco pi- vs pt",100,0,25)
hNum=TH1F("hNum","reco K/pi vs pt",100,0,25)
hDen=TH1F("hDen","gen K/pi vs pt",100,0,25)
hDenK=TH1F("hDenK","gen K vs pt",100,0,25)
hDenKp=TH1F("hDenKp","gen K+ vs pt",100,0,25)
hDenKm=TH1F("hDenKm","gen K- vs pt",100,0,25)
hDenpi=TH1F("hDenpi","gen pi vs pt",100,0,25)
hDenpip=TH1F("hDenpip","gen pi+ vs pt",100,0,25)
hDenpim=TH1F("hDenpim","gen pi- vs pt",100,0,25)
hdR=TH1F("hdR","dR distribution",1000,0,1)            
hdRK=TH1F("hdRK","dR distribution (K)",100,0,0.02)
hdRpi=TH1F("hdRpi","dR distribution (pi)",100,0,1)
hNumPhi=TH1F("hNumPhi","reco Phi vs pt",100,0,25)
hDenPhi=TH1F("hDenPhi","Gen Phi vs pt",100,0,25)
hNumKs=TH1F("hNumKs","reco Ks vs pt",100,0,25)
hDenKs=TH1F("hDenKs","Gen Ks vs pt",100,0,25)
hdRPhi=TH1F("hdRPhi","dR distribution (Phi)",100,0,25)
hdRKs=TH1F("hdRKs","dR distribution (Ks)",100,0,25)
hNumB=TH1F("hNumB","reco B vs pt",100,0,25)
hDenB=TH1F("hDenB","gen B vs pt",100,0,25)
hdRB=TH1F("hdRB","dR distribution (B)",100,0,25)


# In[350]:


def matching1(tree,p,m,hDen,hNum,hdR,th,hIM):
    matches=[]
    for plus in p:
        for minus in m:
            if plus[1]==minus[1]:
                pt=tree.GenPart_pt[plus[1]]
                eta=tree.GenPart_pt[plus[1]]
                phi=tree.GenPart_pt[plus[1]]
                mass=tree.GenPart_mass[plus[1]]
                GenMother=TLorentzVector()
                GenMother.SetPtEtaPhiM(pt,eta,phi,mass)
                RecoMother=plus[2]+minus[2]
                dr=GenMother.DeltaR(RecoMother)
                hIM.Fill(RecoMother.M())
                hdR.Fill(dr)
                if dr<th:
                    hNum.Fill(pt)
                    matches.append((GenMother,tree.GenPart_genPartIdxMother[plus[1]],RecoMother))
    return matches          


# In[351]:


def matching2 (Gp,Tp):
    pdg = ROOT.TDatabasePDG.Instance()
    mK=[]
    mpi=[]
    go=True
    while go:
        go=False
        match=dict()
        i=0
        while i<(len(Gp)):
            dR0=1000
            idx=0
            for j in range(len(Tp)):
                dR=Gp[i][0].DeltaR(Tp[j])
                if dR<dR0:
                    dR0=dR
                    idx=j
            if abs(Gp[i][2])==211 and dR0<0.1:
                if idx in match:
                    go=True
                    match[idx].add((i,dR0))
                else:
                    match[idx]={(i,dR0)}
                i+=1
            elif abs(Gp[i][2])==321 and dR0<0.02:
                if idx in match:
                    go=True
                    match[idx].add((i,dR0))
                else:
                    match[idx]={(i,dR0)}
                i+=1
            else:
                hdR.Fill(dR0)
                if abs(Gp[i][2])==321:
                    hdRK.Fill(dR0)
                if abs(Gp[i][2])==211:
                    hdRpi.Fill(dR0)
                del Gp[i]
        l_g=[]
        l_t=[]
        for e in match:
            dr=1000
            ind=0
            for el in match[e]:
                if el[1]<dr:
                    dr=el[1]
                    ind=el[0]
            hdR.Fill(dr)
            Tp[e].SetPtEtaPhiM(Tp[e].Pt(),Tp[e].Eta(),Tp[e].Phi(),pdg.GetParticle(Gp[ind][2]).Mass())
            if Gp[ind][2]==321:
                hNum.Fill(Gp[ind][3])
                hNumK.Fill(Gp[ind][3])
                hNumKp.Fill(Gp[ind][3])
                hdRK.Fill(dr)
                mK.append((Gp[ind][0],Gp[ind][1],Tp[e]))
            elif Gp[ind][2]==-321:
                hNum.Fill(Gp[ind][3])
                hNumK.Fill(Gp[ind][3])
                hNumKm.Fill(Gp[ind][3])
                hdRK.Fill(dr)
                mK.append((Gp[ind][0],Gp[ind][1],Tp[e]))
            elif Gp[ind][2]==211:
                hNum.Fill(Gp[ind][3])
                hNumpi.Fill(Gp[ind][3])
                hNumpip.Fill(Gp[ind][3])
                hdRpi.Fill(dr)
                mpi.append((Gp[ind][0],Gp[ind][1],Tp[e]))
            else:
                hNum.Fill(Gp[ind][3])
                hNumpi.Fill(Gp[ind][3])
                hNumpim.Fill(Gp[ind][3])
                hdRpi.Fill(dr)
                mpi.append((Gp[ind][0],Gp[ind][1],Tp[e]))
            l_g.append(ind)
            l_t.append(e)
        l_g.sort(reverse=True)
        l_t.sort(reverse=True)
        for e in l_g:
            del Gp[e]
        for e in l_t:
            del Tp[e]
    return (mK,mpi)


# In[352]:


h_IM_Phi.Reset()
h_IM_Ks.Reset()
h_IM_B.Reset()
h_pt_Kp.Reset()
h_pt_Km.Reset()
h_pt_pip.Reset()
h_pt_pim.Reset()
h_Eta_Kp.Reset()
h_Eta_Km.Reset()
h_Eta_pip.Reset()
h_Eta_pim.Reset()
h_Phi_Kp.Reset()
h_Phi_Km.Reset()
h_Phi_pip.Reset()
h_Phi_pim.Reset()
hNumK.Reset()
hNumKp.Reset()
hNumKm.Reset()
hNumpi.Reset()
hNumpip.Reset()
hNumpim.Reset()
hNum.Reset()
hDen.Reset()
hDenK.Reset()
hDenKp.Reset()
hDenKm.Reset()
hDenpi.Reset()
hDenpip.Reset()
hDenpim.Reset()
hdR.Reset()
hdRK.Reset()
hdRpi.Reset()
hNumPhi.Reset()
hDenPhi.Reset()
hNumKs.Reset()
hDenKs.Reset()
hNumB.Reset()
hDenB.Reset()
hdRPhi.Reset()
hdRB.Reset()
hdRKs.Reset()
pdg = ROOT.TDatabasePDG.Instance()
for i in range(N):
    tree.GetEntry(i)
    Kp=[]
    Km=[]
    pip=[]
    pim=[]
    for id,midx,pt,eta,phi,st,mass in zip(tree.GenPart_pdgId,tree.GenPart_genPartIdxMother,tree.GenPart_pt,tree.GenPart_eta,tree.GenPart_phi,tree.GenPart_status,tree.GenPart_mass):
        gmidx=tree.GenPart_genPartIdxMother[midx]
        if abs(tree.GenPart_pdgId[gmidx])==511 and st==1:
            if abs((tree.GenPart_pdgId[midx]))==333:
                if id==321:
                    hDenB.Fill(tree.GenPart_pt[gmidx])
                    hDenPhi.Fill(tree.GenPart_pt[midx])
                    hDen.Fill(pt)
                    hDenK.Fill(pt)
                    hDenKp.Fill(pt)
                    GenP=TLorentzVector()
                    GenP.SetPtEtaPhiM(pt,eta,phi,mass)
                    Kp.append((GenP,midx,id,pt))
                elif id==-321:
                    hDen.Fill(pt)
                    hDenK.Fill(pt)
                    hDenKm.Fill(pt)
                    GenP=TLorentzVector()
                    GenP.SetPtEtaPhiM(pt,eta,phi,mass)
                    Km.append((GenP,midx,id,pt))
            elif abs((tree.GenPart_pdgId[midx]))==310:
                if id==211:
                    hDenKs.Fill(tree.GenPart_pt[midx])
                    hDen.Fill(pt)
                    hDenpi.Fill(pt)
                    hDenpip.Fill(pt)
                    GenP=TLorentzVector()
                    GenP.SetPtEtaPhiM(pt,eta,phi,pdg.GetParticle(id).Mass())
                    pip.append((GenP,midx,id,pt))
                elif id==-211:
                    hDen.Fill(pt)
                    hDenpi.Fill(pt)
                    hDenpim.Fill(pt)
                    GenP=TLorentzVector()
                    GenP.SetPtEtaPhiM(pt,eta,phi,pdg.GetParticle(id).Mass())
                    pim.append((GenP,midx,id,pt))
    Tp=[]
    Tm=[]
    for ch,t_pt,t_eta,t_phi,t_mass in zip(tree.Track_charge,tree.Track_pt,tree.Track_eta,tree.Track_phi,tree.Track_mass):
        Track=TLorentzVector()
        Track.SetPtEtaPhiM(t_pt,t_eta,t_phi,t_mass)
        if ch==1:
            Tp.append(Track)
        elif ch==-1:
            Tm.append(Track)
    mKp,mpip=matching2(Kp+pip,Tp)
    mKm,mpim=matching2(Km+pim,Tm)
    for e in mKp:
        h_pt_Kp.Fill(e[2].Pt())
        h_Eta_Kp.Fill(e[2].Eta())
        h_Phi_Kp.Fill(e[2].Phi())
    for e in mKm:
        h_pt_Km.Fill(e[2].Pt())
        h_Eta_Km.Fill(e[2].Eta())
        h_Phi_Km.Fill(e[2].Phi())
    for e in mpip:
        h_pt_pip.Fill(e[2].Pt())
        h_Eta_pip.Fill(e[2].Eta())
        h_Phi_pip.Fill(e[2].Phi())
    for e in mpim:
        h_pt_pim.Fill(e[2].Pt())
        h_Eta_pim.Fill(e[2].Eta())
        h_Phi_pim.Fill(e[2].Phi())
    mPhi=matching1(tree,mKp,mKm,hDenPhi,hNumPhi,hdRPhi,20,h_IM_Phi)
    mKs=matching1(tree,mpip,mpim,hDenKs,hNumKs,hdRKs,20,h_IM_Ks)
    mB=matching1(tree,mPhi,mKs,hDenB,hNumB,hdRB,20,h_IM_B)


# In[353]:


os.system("mkdir -p Plots")


# In[354]:


c0=TCanvas("c0","pt distributions",200,10,700,500)
c0.cd()
h_pt_Kp.SetLineColor(2)
h_pt_Km.SetLineColor(4)
h_pt_pip.SetLineColor(6)
h_pt_pim.SetLineColor(8)
h_pt_Kp.SetStats(0) 
h_pt_Kp.Draw()
h_pt_Km.Draw("same")
h_pt_pip.Draw("same")
h_pt_pim.Draw("same")
h_pt_Kp.SetTitle("Pt distributions")
h_pt_Kp.GetXaxis().SetTitle("p_{t} [Gev/c]")
h_pt_Kp.GetYaxis().SetTitle("number of tracks")
legend=TLegend(0.7,0.8,0.85,0.6)
legend.AddEntry(h_pt_Kp,"K+ p_{t}")
legend.AddEntry(h_pt_Km,"K- p_{t}")
legend.AddEntry(h_pt_pip,"pi+ p_{t}")
legend.AddEntry(h_pt_pim,"pi- p_{t}")
ROOT.gStyle.SetLegendTextSize(0.05)
legend.Draw()
c0.Draw()
c0.SaveAs("Plots/pt_distributions.png")


# In[355]:


c00=TCanvas("c00","Eta distributions",200,10,700,500)
c00.cd()
h_Eta_Kp.SetLineColor(2)
h_Eta_Km.SetLineColor(4)
h_Eta_pip.SetLineColor(6)
h_Eta_pim.SetLineColor(8)
h_Eta_Kp.SetStats(0)
h_Eta_Kp.Draw()
h_Eta_Km.Draw("same")
h_Eta_pip.Draw("same")
h_Eta_pim.Draw("same")
h_Eta_Kp.SetTitle("Eta distributions")
h_Eta_Kp.GetXaxis().SetTitle("Eta")
h_Eta_Kp.GetYaxis().SetTitle("number of tracks")
legend=TLegend(0.7,0.8,0.85,0.6)
legend.AddEntry(h_Eta_Kp,"K+ Eta ")
legend.AddEntry(h_Eta_Km,"K- Eta")
legend.AddEntry(h_Eta_pip,"pi+ Eta ")
legend.AddEntry(h_Eta_pim,"pi- Eta")
ROOT.gStyle.SetLegendTextSize(0.05)
legend.Draw()
c00.Draw()
c00.SaveAs("Plots/Eta_distributions.png")


# In[356]:


c000=TCanvas("c000","Phi distributions",200,10,700,500)
c000.cd()
h_Phi_Kp.SetLineColor(2)
h_Phi_Km.SetLineColor(4)
h_Phi_pip.SetLineColor(6)
h_Phi_pim.SetLineColor(8)
h_Phi_Kp.SetStats(0)
h_Phi_Kp.Draw()
h_Phi_Km.Draw("same")
h_Phi_pip.Draw("same")
h_Phi_pim.Draw("same")
h_Phi_Kp.SetTitle("Eta distributions")
h_Phi_Kp.GetXaxis().SetTitle("Eta")
h_Phi_Kp.GetYaxis().SetTitle("number of tracks")
legend=TLegend(0.7,0.8,0.85,0.6)
legend.AddEntry(h_Phi_Kp,"K+ Phi ")
legend.AddEntry(h_Phi_Km,"K- Phi")
legend.AddEntry(h_Phi_pip,"pi+ Phi ")
legend.AddEntry(h_Phi_pim,"pi- Phi")
ROOT.gStyle.SetLegendTextSize(0.05)
legend.Draw()
c000.Draw()
c000.SaveAs("Plots/Phi_distributions.png")


# In[357]:


c0000=TCanvas("c0000","Phi invariant mass",200,10,700,500)
c0000.cd()
h_IM_Phi.Draw()
h_IM_Phi.GetXaxis().SetTitle("invariant mass [Gev/c]")
h_IM_Phi.GetYaxis().SetTitle("number of couples K+K-")
c0000.Draw()
c0000.SaveAs("Plots/K+K- invariant mass.png")


# In[358]:


c00000=TCanvas("c00000","Ks invariant mass",200,10,700,500)
c00000.cd()
h_IM_Ks.Draw()
h_IM_Ks.GetXaxis().SetTitle("invariant mass [Gev/c]")
h_IM_Ks.GetYaxis().SetTitle("number of couples pi+pi-")
c00000.Draw()
c00000.SaveAs("Plots/pi+pi- invariant mass.png")


# In[359]:


c000000=TCanvas("c000000","Phi invariant mass",200,10,700,500)
c000000.cd()
h_IM_B.Draw()
h_IM_B.GetXaxis().SetTitle("invariant mass [Gev/c]")
h_IM_B.GetYaxis().SetTitle("number of quartets K+K-pi+pi-")
c000000.Draw()
c000000.SaveAs("Plots/K+K-pi+pi- invariant mass.png")


# In[360]:


h=TH1F("h","Efficiency(K) vs pt",100,0,25)
h.Reset()
c1=TCanvas('c1','Efficiency K vs Pt',200,10,700,500)
c1.cd()
EffK_pt=TEfficiency(hNumK,hDenK)
h.GetXaxis().SetTitle("Gen K p_{t} [Gev/c]")
h.GetYaxis().SetTitle("Efficiency")
ROOT.gStyle.SetOptStat("e")
h.SetStats(0)
h.Draw()
EffK_pt.Draw("same")
c1.Draw()
c1.SaveAs("Plots/Efficiency(K)_vs_pt.png")


# In[361]:


h1=TH1F("h1","Efficiency(K+) vs pt",100,0,25)
h1.Reset()
c1a=TCanvas('c1a','Efficiency K+ vs Pt',200,10,700,500)
c1a.cd()
EffKp_pt=TEfficiency(hNumKp,hDenKp)
h1.GetXaxis().SetTitle("Gen K+ p_{t} [Gev/c]")
h1.GetYaxis().SetTitle("Efficiency")
h1.SetStats(0)
h1.Draw()
EffKp_pt.Draw("same")
c1a.Draw()
c1a.SaveAs("Plots/Efficiency(K+)_vs_pt.png")


# In[362]:


h1b=TH1F("h1b","Efficiency(K-) vs pt",100,0,25)
h1b.Reset()
c1b=TCanvas('c1b','Efficiency K- vs Pt',200,10,700,500)
c1b.cd()
EffKm_pt=TEfficiency(hNumKm,hDenKm)
h1b.GetXaxis().SetTitle("Gen K- p_{t} [Gev/c]")
h1b.GetYaxis().SetTitle("Efficiency")
h1b.SetStats(0)
h1b.Draw()
EffKm_pt.Draw("same")
c1b.Draw()
c1b.SaveAs("Plots/Efficiency(K-)_vs_pt.png")


# In[363]:


h2=TH1F("h2","Efficiency(pi) vs pt",100,0,25)
h2.Reset()
c2=TCanvas('c2','Efficiency pi vs Pt',200,10,700,500)
c2.cd()
Effpi_pt=TEfficiency(hNumpi,hDenpi)
h2.GetXaxis().SetTitle("Gen pi p_{t} [Gev/c]")
h2.GetYaxis().SetTitle("Efficiency")
h2.SetStats(0)
h2.Draw()
Effpi_pt.Draw("same")
c2.Draw()
c2.SaveAs("Plots/Efficiency(pi)_vs_pt.png")


# In[364]:


h2a=TH1F("h2a","Efficiency(pi+) vs pt",100,0,25)
h2a.Reset()
c2a=TCanvas('c2a','Efficiency pi+ vs Pt',200,10,700,500)
c2a.cd()
Effpip_pt=TEfficiency(hNumpip,hDenpip)
h2a.GetXaxis().SetTitle("Gen pi+ p_{t} [Gev/c]")
h2a.GetYaxis().SetTitle("Efficiency")
h2a.SetStats(0)
h2a.Draw()
Effpip_pt.Draw("same")
c2a.Draw()
c2a.SaveAs("Plots/Efficiency(pi+)_vs_pt.png")


# In[365]:


h2b=TH1F("h2b","Efficiency(pi-) vs pt",100,0,25)
h2b.Reset()
c2b=TCanvas('c2b','Efficiency pi- vs Pt',200,10,700,500)
c2b.cd()
Effpim_pt=TEfficiency(hNumpim,hDenpim)
h2b.GetXaxis().SetTitle("Gen pi- p_{t} [Gev/c]")
h2b.GetYaxis().SetTitle("Efficiency")
h2b.SetStats(0)
h2b.Draw()
Effpim_pt.Draw("same")
c2b.Draw()
c2b.SaveAs("Plots/Efficiency(pi-)_vs_pt.png")


# In[366]:


h3=TH1F("h3","Efficiency vs pt",100,0,25)
h3.Reset()
c3=TCanvas('c3','Efficiency_Pt',200,10,700,500)
c3.cd()
Eff_pt=TEfficiency(hNum,hDen)
h3.GetXaxis().SetTitle("Gen pi/K p_{t} [Gev/c]")
h3.GetYaxis().SetTitle("Efficiency")
h3.SetStats(0)
h3.Draw()
Eff_pt.Draw("same")
c3.Draw()
c3.SaveAs("Plots/Efficiency_vs_pt.png")
    


# In[367]:


c5=TCanvas("c5","K",200,10,700,500)
c5.cd()
hNumK.SetLineColor(2)
hDenK.SetLineColor(4)
hNumK.Draw()
hDenK.Draw("same")
hNumK.SetTitle("Gen K compared to matched gen K vs pt")
hNumK.GetXaxis().SetTitle("Gen K p_{t} [Gev/c]")
hNumK.GetYaxis().SetTitle("number of Kaons")
legend=TLegend(0.6,0.8,0.85,0.6)
legend.AddEntry(hNumK,"matched gen K")
legend.AddEntry(hDenK,"gen K")
ROOT.gStyle.SetLegendTextSize(0.05)
legend.Draw()
c5.Draw()
c5.SaveAs("Plots/Gen_vs_Reco_K.png")


# In[368]:


c6=TCanvas("c6","pi",200,10,700,500)
c6.cd()
hNumpi.SetLineColor(2)
hDenpi.SetLineColor(4)
hNumpi.Draw()
hDenpi.Draw("same")
hNumpi.SetTitle("Gen pi compared to matched gen pi vs pt")
hNumpi.GetXaxis().SetTitle("Gen pi p_{t} [Gev/c]")
hNumpi.GetYaxis().SetTitle("number of pions")
legend=TLegend(0.6,0.8,0.85,0.6)
legend.AddEntry(hNumpi,"matched gen pi")
legend.AddEntry(hDenpi,"gen pi")
ROOT.gStyle.SetLegendTextSize(0.05)
legend.Draw()
c6.Draw()
c6.SaveAs("Plots/Gen_vs_Reco_pi.png")


# In[369]:


c9=TCanvas("c9","dR",200,10,700,500)
hdR.Draw()
c9.Draw()
c9.SaveAs("Plots/dR.png")


# In[370]:


c10=TCanvas("c10","dR dei K",200,10,700,500)
hdRK.Draw()
c10.Draw()
c10.SaveAs("Plots/dR(K).png")


# In[371]:


c11=TCanvas("c11","dR dei pi",200,10,700,500)
hdRpi.Draw()
c11.Draw()
c11.SaveAs("Plots/dR(pi).png")


# In[372]:


h10=TH1F("h10","Efficiency(Phi) vs pt",100,0,25)
h10.Reset()
c12=TCanvas('c12','Efficiency Phi vs Pt',200,10,700,500)
c12.cd()
EffPhi_pt=TEfficiency(hNumPhi,hDenPhi)
h10.GetXaxis().SetTitle("Gen Phi p_{t} [Gev/c]")
h10.GetYaxis().SetTitle("Efficiency")
h10.SetStats(0)
h10.Draw()
EffPhi_pt.Draw("same")
c12.Draw()
c12.SaveAs("Plots/Efficiency(Phi)_vs_pt.png")


# In[373]:


h11=TH1F("h11","Efficiency(Ks) vs pt",100,0,25)
h11.Reset()
c13=TCanvas('c13','Efficiency Ks vs Pt',200,10,700,500)
c13.cd()
EffKs_pt=TEfficiency(hNumKs,hDenKs)
h11.GetXaxis().SetTitle("Gen Phi p_{t} [Gev/c]")
h11.GetYaxis().SetTitle("Efficiency")
h11.SetStats(0)
h11.Draw()
EffKs_pt.Draw("same")
c13.Draw()
c13.SaveAs("Plots/Efficiency(Ks)_vs_pt.png")


# In[374]:


c14=TCanvas("c14","dR dei Phi",200,10,700,500)
hdRPhi.Draw()
c14.Draw()
c14.SaveAs("Plots/dR(Phi).png")


# In[375]:


c15=TCanvas("c15","dR dei Ks",200,10,700,500)
hdRKs.Draw()
c15.Draw()
c15.SaveAs("Plots/dR(Ks).png")


# In[376]:


h12=TH1F("h12","Efficiency(B) vs pt",100,0,25)
h12.Reset()
c16=TCanvas('c16','Efficiency B vs Pt',200,10,700,500)
c16.cd()
EffB_pt=TEfficiency(hNumB,hDenB)
h12.GetXaxis().SetTitle("Gen Phi p_{t} [Gev/c]")
h12.GetYaxis().SetTitle("Efficiency")
h12.SetStats(0)
h12.Draw()
EffB_pt.Draw("same")
c16.Draw()
c16.SaveAs("Plots/Efficiency(B)_vs_pt.png")


# In[377]:


c17=TCanvas("c17","dR dei B",200,10,700,500)
hdRB.Draw()
c17.Draw()
c17.SaveAs("Plots/dR(B).png")


# In[378]:


get_ipython().system('jupyter nbconvert --to script Untitled.ipynb')


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




