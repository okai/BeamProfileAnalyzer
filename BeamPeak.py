#!/usr/bin/env python
# -*- coding: utf-8 -*-


import tkinter
import tkinter.filedialog as tkfd
import h5py
import numpy as np
import matplotlib
import sys
import time
from mpl_toolkits.mplot3d import Axes3D
import array
import struct
import io
import subprocess
import threading
import datetime
import os
import os.path
matplotlib.use('TkAgg')
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt


class ReaderApp( tkinter.Frame ):
    def __init__( self, master=None ):
        tkinter.Frame.__init__( self, master )
        self.grid()
        self.createConstants()
        self.createWidgets()

    def createConstants( self ):
        self.power = tkinter.StringVar()
        self.smoothing = tkinter.StringVar()

    def createWidgets( self ):

        pathToDataButton = tkinter.Button( self, text=u"   choose bgData   " , command=self.AskFileName )
        pathToDataButton.grid( row=1, column=0 )

        powerLabel = tkinter.Label( self, text=u"power " )
        powerEntry = tkinter.Entry( self, textvariable=self.power )
        powerUnitLabel = tkinter.Label( self, text=u" mJ")
        powerLabel.grid( row=2, column=0 )
        powerEntry.grid( row=2, column=1 )
        powerUnitLabel.grid( row=2, column=2 )

        smoothingLabel = tkinter.Label( self, text=u"smoothing" )
        smoothingEntry = tkinter.Entry( self, textvariable=self.smoothing )
        smoothingLabel.grid( row=3, column=0 )
        smoothingEntry.grid( row=3, column=1, columnspan=2 )

        startButton = tkinter.Button( self, text="OK", command=self.Acquisition )
        startButton.grid( row=4, column=0 )

    def AskFileName( self ):
        self.pathToData = tkfd.askopenfilename()
        print(self.pathToData)
        pathToDataLabel = tkinter.Label( self, text=os.path.basename(self.pathToData))
        pathToDataLabel.grid( row=1, column=1 )


    def Acquisition( self ):

        t1 = time.time()
        #########################################################################
        # bgData
        FilePath = self.pathToData
        FileName = os.path.basename(FilePath)

        DatName = os.path.splitext(FileName)[0]
        # Beam power [mJ]
        Pwr = float(self.power.get())

        # Mesh size for Smoothing
        S = int(self.smoothing.get())

        # Output directory
        Dir = os.path.splitext(FilePath)[0]

        #########################################################################
        print("\nDirectory =", Dir)
        print('\nFile = %s.bgData' % (DatName))
        print('\nBeam power = %.1f [mJ]' % (Pwr))

        #########################################################################
        # output directory
        if not os.path.exists(Dir):
            os.mkdir(Dir)

        # output file names
        Out1 = Dir + '/Data_array.txt'
        Out2 = Dir + '/Smoothed_Data_array.txt'
        Out3 = Dir + '/Projection_Y.txt'
        Out4 = Dir + '/Projection_X.txt'
        Out5 = Dir + '/Section_X_at_Max.txt'
        Out6 = Dir + '/Section_Y_at_Max.txt'

        #########################################################################
        # the procedure describing blow is to read the contents of file name 1
        BgData = h5py.File(FilePath, 'r')

        H = BgData['BG_DATA/1/RAWFRAME/HEIGHT'].value
        W = BgData['BG_DATA/1/RAWFRAME/WIDTH'].value

        BinX = BgData['BG_DATA/1/RAWFRAME/BINNINGX'].value
        BinY = BgData['BG_DATA/1/RAWFRAME/BINNINGY'].value
        PX = BinX * BgData['BG_DATA/1/RAWFRAME/PIXELSCALEXUM'].value  # μm
        PY = BinY * BgData['BG_DATA/1/RAWFRAME/PIXELSCALEYUM'].value  # μm

        print(H[0], W[0], BinX[0], BinY[0], PX[0], PY[0])

        GS = BgData['BG_DATA/1//RAWFRAME/GAINSTAMP'].value
        BE = BgData['BG_DATA/1/RAWFRAME/BITENCODING'].value
        BeD = BE[0].decode('utf-8')

        if BeD == 'L8':
            CF = (2 ** (32 - 1 - 8))
        elif BeD == 'L16_8':
            CF = (2 ** (32 - 1 - 8))
        elif BeD == 'L16_10':
            CF = (2 ** (32 - 1 - 10))
        elif BeD == 'L16_12':
            CF = (2 ** (32 - 1 - 12))
        else:
            sys.exit()

        Raw = BgData['BG_DATA/1/DATA'].value
        Cal = Raw / CF

        SO = BgData['BG_DATA/1/RAWFRAME/SUBTRACTION_OFFSET'].value

        ################################################################################
        ##Converting the date of file1, file2 and subtraction results to the H*W matrix
        Sum = np.sum(Cal)
        DatR = Pwr * Cal / Sum / PX[0] / PY[0] * 10 ** 5
        Dat = DatR.reshape([H[0], W[0]])

        WR = PX[0] * int(W[0]) / 1000  # mm
        HR = PY[0] * int(H[0]) / 1000  # mm
        print('FrameSize = %.2f mm * %.2f mm' % (WR, HR))

        #############################################################################
        ##Projection of Dat onto x and y >>> YPr, XPr
        YPr = np.sum(Dat, axis=1)
        XPr = np.sum(Dat, axis=0)

        Ny = np.arange(H[0]) * PY[0] / 1000
        Nx = np.arange(W[0]) * PX[0] / 1000

        # YPr_dev = np.std(YPr)
        # XPr_dev = np.std(XPr)
        ################################################################################
        # Smoothing Dat >>> SDat
        H2 = int(H[0] / S)
        W2 = int(W[0] / S)

        SDat = np.zeros((H2, W2))

        for i in range(H2):
            for j in range(W2):
                SDat[i, j] = np.sum(Dat[(S * i):(S * i + S), (S * j):(S * j + S)]) / (S ** 2)

        SNy = np.arange(H2) * PY[0] * S / 1000
        SNx = np.arange(W2) * PX[0] * S / 1000

        ################################################################################
        # Obteining the position of peak
        Max = np.max(DatR)
        MaxP = np.where(Dat == Max)
        MaxPx = MaxP[1] * PX[0] / 1000
        MaxPy = MaxP[0] * PY[0] / 1000

        SMax = SDat.max()
        SMaxP = np.where(SDat == SMax)
        SMaxPx = SMaxP[1] * PX[0] * S / 1000
        SMaxPy = SMaxP[0] * PY[0] * S / 1000

        # print('\nMax(X,Y) = ( %.3f , %.3f ) = %.3f '%(MaxPx, MaxPy, Max))
        print('Density =  %.3f [J/cm2]' % (Max))

        print('\nSmoothed = %d * %d ' % (S, S))
        # print('\nMax(X,Y) = ( %.3f , %.3f ) = %.3f '%(SMaxPx, SMaxPy, SMax))
        print('Density =  %.3f [J/cm2]' % (SMax))

        ################################################################################
        # Getting segments of smoothed data
        XSc = SDat[SMaxP[0], :]
        YSc = SDat[:, SMaxP[1]]

        ##############################################################################
        ##Saving the result to the txt files
        np.savetxt(Out1, Dat, fmt='%0.2f')
        np.savetxt(Out3, YPr.T, fmt='%0.2f')
        np.savetxt(Out4, XPr.T, fmt='%0.2f')

        np.savetxt(Out2, SDat, fmt='%0.2f')
        np.savetxt(Out5, XSc.T, fmt='%0.2f')
        np.savetxt(Out6, YSc, fmt='%0.2f')

        ##############################################################################
        # Prepairing field for figures.
        fig = plt.figure(figsize=(14, 8.3))
        gs = gridspec.GridSpec(4, 4, width_ratios=[4, 1, 4, 1], height_ratios=[5, 15, 1, 3])

        axRaw = plt.subplot(gs[1, 0])
        axProX = plt.subplot(gs[0, 0], sharex=axRaw)
        axProY = plt.subplot(gs[1, 1], sharey=axRaw)

        axS = plt.subplot(gs[1, 2])
        axScX = plt.subplot(gs[0, 2], sharex=axS)
        axScY = plt.subplot(gs[1, 3], sharey=axS)

        axCB_1 = plt.subplot(gs[2, 0])
        axCB_2 = plt.subplot(gs[2, 2])

        axTxt_1 = plt.subplot(gs[3, 0:2])
        axTxt_2 = plt.subplot(gs[3, 2:4])

        ###############################################################################
        ##Displaying the distribution of intensity
        Fig1 = axRaw.imshow(Dat, interpolation='none', origin='lower', extent=[0, WR, 0, HR])
        Fig1.set_cmap('nipy_spectral')
        cb_1 = plt.colorbar(Fig1, orientation='horizontal', cax=axCB_1)
        cb_1.set_label('Intensity [J/mm2]')
        axRaw.set_title(DatName)
        axRaw.set_xlabel('X [mm]')
        axRaw.set_ylabel('Y [mm]')
        plt.tight_layout()

        ###############################################################################
        ##Displaying the image of smoothed distribution of intensity
        Fig2 = axS.imshow(SDat, interpolation='none', origin='lower', extent=[0, WR, 0, HR])
        Fig2.set_cmap('nipy_spectral')
        cb_2 = plt.colorbar(Fig2, orientation='horizontal', cax=axCB_2)
        cb_2.set_label('Intensity [J/mm2]')
        axS.set_title('Smoothed')
        axS.set_xlabel('X [mm]')
        axS.set_ylabel('Y [mm]')
        plt.tight_layout()

        ###############################################################################
        ##Ploting the x y projection
        Fig4 = axProY.plot(YPr, Ny)
        axProY.ticklabel_format(style='sci', scilimits=(0, 0), axis='x')
        axProY.set_title('Y-Projection')
        axProY.set_ylabel('Y [mm]')
        axProY.set_xlabel('Count')
        plt.tight_layout()

        Fig5 = axProX.plot(Nx, XPr)
        axProX.ticklabel_format(style='sci', scilimits=(0, 0), axis='y')
        axProX.set_title('X-projection')
        axProX.set_xlabel('X [mm]')
        axProX.set_ylabel('Count')
        plt.tight_layout()

        ###############################################################################
        ##Ploting the x y section
        Fig6 = axScY.plot(YSc.T[0, :], SNy)
        axScY.ticklabel_format(style='sci', scilimits=(0, 0), axis='x')
        axScY.set_title('Y-Section')
        axScY.set_ylabel('Y [mm]')
        axScY.set_xlabel('Intensity \n [J/cm2]')
        plt.tight_layout()

        Fig7 = axScX.plot(SNx, XSc[0, :])
        axScX.ticklabel_format(style='sci', scilimits=(0, 0), axis='y')
        axScX.set_title('X-Section')
        axScX.set_xlabel('X [mm]')
        axScX.set_ylabel('Intensity \n [J/cm2]')
        plt.tight_layout()

        ###############################################################################
        ##text
        # Txt_1 = '\nBeam power = %.1f mJ\n' % (Pwr) \
        #         + 'Max(X,Y) = ( %.3f mm , %.3f mm ) = %.3f \n' % (MaxPx, MaxPy, Max) \
        #         + 'Density = %.3f J / cm2' % (Max)
        #
        # Txt_2 = 'Frame size = %.3f mm * %.3f mm \n' % (WR, HR) \
        #         + 'Smoothed = %d * %d \n' % (S, S) \
        #         + 'Max(X,Y) = ( %.3f mm , %.3f mm ) = %.3f \n' % (SMaxPx, SMaxPy, SMax) \
        #         + 'Density = %.3f J / cm2' % (SMax)
        #
        # axTxt_1.axis('off')
        # axTxt_2.axis('off')
        # axTxt_1.text(0, 0.2, Txt_1, fontsize=11)
        # axTxt_2.text(0, 0.2, Txt_2, fontsize=11)

        ###############################################################################
        ##save
        savePDF = Dir + '/' + DatName + '.pdf'
        savePNG = Dir + '/' + DatName + '.png'

        fig.savefig(savePDF)
        fig.savefig(savePNG)

        ###############################################################################
        t2 = time.time()
        elapsed_time = t2 - t1
        print('\n time : %f sec' % (elapsed_time))

        ###############################################################################
        plt.show()

        print('大丈夫ですか？')

if __name__=='__main__':
    app = ReaderApp()
    app.master.title( "おつかれさまです！！！")
    app.mainloop()

# app.SDD_raw.close()

