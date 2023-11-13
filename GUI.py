# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'Data3.ui'
#
# Created by: PyQt5 UI code generator 5.9.2
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import QFileDialog
from PyQt5.QtWidgets import QMessageBox
from PyQt5.QtWidgets import QLineEdit
import sys
from PyQt5.QtWidgets import QApplication, QWidget, QInputDialog, QComboBox
from PyQt5.QtGui import QIcon
import os
import csv
import numpy as np
from scipy import stats
import pandas as pd
import matplotlib.pyplot as plt
import shutil

class Ui_MainWindow(object):
    
    xmin=4000
    xmax=500
    legend=2   
        
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(800, 600)
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.label = QtWidgets.QLabel(self.centralwidget)
        self.label.setGeometry(QtCore.QRect(200, 20, 371, 61))
        font = QtGui.QFont()
        font.setPointSize(18)
        font.setBold(False)
        font.setWeight(50)
        self.label.setFont(font)
        self.label.setObjectName("label")
        
        self.ff = 0
        self.radioButton_file = QtWidgets.QRadioButton(self.centralwidget)
        self.radioButton_file.setGeometry(QtCore.QRect(70, 100, 95, 20))
        self.radioButton_file.toggled.connect(self.fileselected)
        
        
        # self.combo = QComboBox()
        # self.combo.setGeometry(QtCore.QRect(70, 100, 371, 41))
        # self.combo.setObjectName(("combo"))
        # self.combo.addItems(['Select File', 'Select Folder'])
        
        self.textBrowser = QtWidgets.QTextBrowser(self.centralwidget)
        self.textBrowser.setGeometry(QtCore.QRect(70, 150, 371, 41))
        self.textBrowser.setObjectName("textBrowser")
        
        self.Browsebutton = QtWidgets.QPushButton(self.centralwidget)
        self.Browsebutton.setGeometry(QtCore.QRect(490, 150, 151, 41))
        self.Browsebutton.setObjectName("Browsebutton")
        
        self.checkBox = QtWidgets.QCheckBox(self.centralwidget)
        self.checkBox.setGeometry(QtCore.QRect(80, 360, 261, 20))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.checkBox.setFont(font)
        
        self.checkBox.setObjectName("checkBox")
        self.checkBox_2 = QtWidgets.QCheckBox(self.centralwidget)
        self.checkBox_2.setGeometry(QtCore.QRect(80, 390, 161, 20))
        font = QtGui.QFont()
        font.setPointSize(10)      
        self.checkBox_2.setFont(font)
        
        self.checkBox_2.setObjectName("checkBox_2")
        self.checkBox_3 = QtWidgets.QCheckBox(self.centralwidget)
        self.checkBox_3.setGeometry(QtCore.QRect(80, 420, 281, 20))
        font = QtGui.QFont()
        font.setPointSize(10)  
        self.checkBox_3.setFont(font)
        
        self.checkBox_3.setObjectName("checkBox_3")
        self.checkBox_4 = QtWidgets.QCheckBox(self.centralwidget)
        self.checkBox_4.setGeometry(QtCore.QRect(80, 450, 181, 20))
        font = QtGui.QFont()
        font.setPointSize(10)
        
        
        self.checkBox_4.setFont(font)
        self.checkBox_4.setObjectName("checkBox_4")
        self.submitbutton = QtWidgets.QPushButton(self.centralwidget)
        self.submitbutton.setGeometry(QtCore.QRect(220, 240, 141, 51))
        
        self.submitbutton.setObjectName("submitbutton")
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 800, 26))       
        
        self.lineEdit_xmin = QtWidgets.QLineEdit(self.centralwidget)
        self.lineEdit_xmin.setGeometry(QtCore.QRect(710, 260, 113, 22))
        self.lineEdit_xmin.setObjectName("lineEdit_xmin")
        self.lineEdit_xmax = QtWidgets.QLineEdit(self.centralwidget)
        self.lineEdit_xmax.setGeometry(QtCore.QRect(710, 320, 113, 22))
        self.lineEdit_xmax.setObjectName("lineEdit_xmax")
        self.xmin_value = QtWidgets.QLabel(self.centralwidget)
        self.xmin_value.setGeometry(QtCore.QRect(710, 230, 91, 16))
        self.xmin_value.setObjectName("xmin_value")
        self.Xmax_value = QtWidgets.QLabel(self.centralwidget)
        self.Xmax_value.setGeometry(QtCore.QRect(710, 290, 91, 16))
        self.Xmax_value.setObjectName("Xmax_value")
        self.Legend = QtWidgets.QLabel(self.centralwidget)
        self.Legend.setGeometry(QtCore.QRect(710, 350, 55, 16))
        self.Legend.setObjectName("Legend")
        self.lineEdit_legend = QtWidgets.QLineEdit(self.centralwidget)
        self.lineEdit_legend.setGeometry(QtCore.QRect(710, 380, 113, 22))
        self.lineEdit_legend.setObjectName("lineEdit_legend")
        self.xlabel = QtWidgets.QLabel(self.centralwidget)
        self.xlabel.setGeometry(QtCore.QRect(710, 410, 91, 16))
        self.xlabel.setObjectName("xlabel")
        self.lineEdit_xlabel = QtWidgets.QLineEdit(self.centralwidget)
        self.lineEdit_xlabel.setGeometry(QtCore.QRect(710, 430, 113, 22))
        self.lineEdit_xlabel.setObjectName("lineEdit_xlabel")
        self.ylabel = QtWidgets.QLabel(self.centralwidget)
        self.ylabel.setGeometry(QtCore.QRect(710, 460, 91, 16))
        self.ylabel.setObjectName("ylabel")
        self.lineEdit_ylabel = QtWidgets.QLineEdit(self.centralwidget)
        self.lineEdit_ylabel.setGeometry(QtCore.QRect(710, 480, 113, 22))
        self.lineEdit_ylabel.setObjectName("lineEdit_ylabel")
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 915, 26))   
        
         
    
        
        self.menubar.setObjectName("menubar")
        MainWindow.setMenuBar(self.menubar)
        
        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)
        

        # Connect Go! push botton to Echem_analysis method
        self.submitbutton.clicked.connect(self.connect_widgets)
        
            
        # Connect file path browse button to get files method
        self.Browsebutton.clicked.connect(self.open_folder_dialog)
    
    def fileselected(self, selected):
        if selected:
            self.ff = 1
            self.label.setText("Select File")
        else:
            self.ff = 0
            self.label.setText("Select Folder")

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "MainWindow"))
        self.label.setText(_translate("MainWindow", "Data Analyser"))
        self.textBrowser.setStatusTip(_translate("MainWindow", "Submit your folder path"))
        self.Browsebutton.setStatusTip(_translate("MainWindow", "Browse folder in PC"))
        self.Browsebutton.setText(_translate("MainWindow", "Browse"))
        self.checkBox.setText(_translate("MainWindow", "Visualize Data in a plot"))
        self.checkBox_2.setText(_translate("MainWindow", "Baseline Corrections"))
        self.checkBox_3.setText(_translate("MainWindow", "Find Absorbance in a specific Wavelength"))
        self.checkBox_4.setText(_translate("MainWindow", "Linear Regreesiion"))
        self.submitbutton.setStatusTip(_translate("MainWindow", "Submit Filepath"))
        self.submitbutton.setText(_translate("MainWindow", "Submit"))
        
        self.xmin_value.setText(_translate("MainWindow", "Enter X min"))
        self.Xmax_value.setText(_translate("MainWindow", "Enter X max"))
        self.Legend.setText(_translate("MainWindow", "Legend"))
        self.xlabel.setText(_translate("MainWindow", "xlabel"))
        self.ylabel.setText(_translate("MainWindow", "ylabel"))
        
        self.radioButton_file.setText(_translate("MainWindow", "Single File"))
        #self.combo.setText(_translate("MainWindow", "Select File"))
    
    
  
    
    
    
    
    
          
    def open_folder_dialog(self):
        if self.ff == 1:
            file_path, _ = QFileDialog.getOpenFileName(caption = 'Single File', 
                                                       directory = os.getcwd(), 
                                                       filter = '*.csv')
            self.textBrowser.setText(file_path)
            
        else:
            
            folder_path = QFileDialog.getExistingDirectory()
            files=QFileDialog.getOpenFileContent()
            print(files)
            self.textBrowser.setText(folder_path)
    
        


  
    
        #print(folder_path)
    # Method that defines folder path, analysis options, and runs Echem_analysis method
    def connect_widgets(self):
    
        # Define file path based on text entered into lineEdit widget
        path = str(self.textBrowser.toPlainText())
    
        # Conditional statement for checking of save plot? check box
        if self.checkBox.isChecked() == True:
            if self.ff == 1:
                self.label.setText("Err: Select Directory!")
            else:
                class_inst = Ui_MainWindow()
                Ui_MainWindow.line_plot(self, path)
                self.label.setText("Data Analysis")
            
        if self.checkBox_2.isChecked() == True:
            if self.ff == 1:
                self.label.setText("Err: Select Directory!")
            else:
                class_inst = Ui_MainWindow()
                Ui_MainWindow.baseline_sub(self, path)
                self.label.setText("Data Analysis")
         
        if self.checkBox_3.isChecked() == True:
            if self.ff == 1:
                self.label.setText("Err: Select Directory!")
            else:
                class_inst = Ui_MainWindow()
                Ui_MainWindow.find_abs(self, path)
                self.label.setText("Data Analysis")
        
            
        if self.checkBox_4.isChecked() == True:
            if self.ff == 0:
                self.label.setText("Err: Select Single!")
            elif self.ff == 1:
                class_inst = Ui_MainWindow()
                Ui_MainWindow.linear_reg(self, path)
                self.label.setText("Data Analysis")
            
        if self.lineEdit_xmin.textChanged() == True:
            Ui_MainWindow.xmin(self, path)
            
    def popup(self):
        msg= QMessageBox()
        msg.setText('Save')
        msg.setWindowTitle('Save?')
        x=msg.exec_()
        
        
    def xmin(self):
        xmin = self.lineEdit_xmin.text()
        print(xmin)
        
    def line_plot(self,path):
        
        xmin=4000
        xmax=500
        legend=0
        xlabel='Wavenumber(cm-1)'
        ylabel='Absorbance'
        #print(path)    
        ### Fetch all files in path
        fileNames = sorted(os.listdir(path))
        
        #print(fileNames)
        
        ### Filter file name list for files ending with .csv
        fileNames = [file for file in fileNames if '.CSV' in file]
        
        ### Predefine empty list to later be filled with max values
        max_value_list = []
        
        ### Loop over all files
        for file in fileNames:
        
            ### Unpack x and y values for each file
            x,y = np.loadtxt(path + '\\' + file, delimiter=',', skiprows = 3, unpack=True)
        
            ### Find max value for each file and append to list
            max_value_list.append(np.max(y))
        
            ### Create line for every file
            plt.plot(x,y,label = file.replace('.CSV',''))
        
        ### x limits
        plt.xlim(xmin,xmax)
        
        ### X-axis label
        plt.xlabel(xlabel, fontsize = 15)
        ### Y-axis label
        plt.ylabel(ylabel, fontsize = 15)
        ### Change tick label (number) size
        plt.tick_params(labelsize = 12)
        
        
        legend=0
        
        ### Options for legend position
        if legend == 0:
            plt.legend(loc='upper right')
        elif legend == 1:
            plt.legend(loc='upper center')
        elif legend == 2:
            plt.legend(loc='upper left')
        else:
            plt.legend(loc = (1.01, 0.35))
        plt.tight_layout()
        
        ### Determine whether to just show plot or save it based on user input
        if saveyn == 'yes':
            plt.savefig(path + '\\' + os.path.basename(path) + '.jpeg', dpi = 300)
        else:
            plt.show()
    
    
    def baseline_sub(self,path):
        
        '''
        Importing multiple csv files, baseline subtracting, plotting as a line, saving as a CSV
        
        PATH = string to path to folder containing multiple files... syntax = r'~file/'
        
        xmin = minimum x value for plotting - integer
        
        xmax = maximum x value for plotting - integer
        
        legend = 0,1, or 2 for legend to be located upper-right, upper-center, or upper-left, respectively
        
        save_plot_yn = 'yes' to save plots, 'no' or anything else to not save plots - string
        
        save_csv_yn = 'yes' to save SNV normalized spectra as CSV, 'no' to not baseline-subtracted
        '''
        
        xmin=4000
        xmax=500
        legend=0 
        ### Fetch all files in path
        fileNames = sorted(os.listdir(path))
        
        ### Filter file name list for files ending with .csv
        fileNames = [file for file in fileNames if '.CSV' in file]
        
        ### Loop over all files
        for file in fileNames:
        
            ### Unpack x and y values for each file
            x,y = np.loadtxt(path + '\\' + file, delimiter=',', skiprows = 3, unpack=True)
        
            ### Find mean at high energy part of reach spectrum (~3800-4000 cm^-1), away from IR resonances
            mean = np.mean(y[885:930])
            
            ### Subtract mean of high energy part of spectrum from entire spectrum, not including spurious points, to subtract baseline
            y_sub = y[73:935] - mean
        
            ### Create line for every unit-normalized spectrum
            plt.plot(x[73:935], y_sub,label = file.replace('.CSV',''))
            
            ### Save csv with each baseline substracted spectrum if save_csv_yn == 'yes'
            #if save_csv_yn == 'yes':
               # x_csv = x[73:935]
               # np.savetxt(path + '\\' + file.replace('.CSV','_sub.CSV'), np.column_stack((x_csv,y_sub)), delimiter = ',')
            #else:
               # print('Type "yes" (note: must be string) for save_csv_yn to save baseline subtracted spectra as their own CSVs')
        
        ### Generate the plot
        plt.xlim(xmin,xmax)
        
        ### X-axis label
        plt.xlabel('Wavenumber (cm^-1)', fontsize = 15)
        ### Y-axis label
        plt.ylabel('Absorbance - baseline subtracted', fontsize = 15)
        ### Change tick label (number) size
        plt.tick_params(labelsize = 12)
        ### Plot title
        plt.title('Baseline subtracted', fontsize = 15)
        ### Call tight_layout so all axis labels are visible in the saved image
        plt.tight_layout()
        
        ### Options for legend position
        if legend == 0:
            plt.legend(loc='upper right')
        elif legend == 1:
            plt.legend(loc='upper center')
        elif legend == 2:
            plt.legend(loc='upper left')
        
        ### Determine whether to just show plot or save it based on user input
        #if save_plot_yn == 'yes':
           # plt.tight_layout()
            #plt.savefig(path + '\\' + os.path.basename(PATH) + '.jpeg', dpi = 300)
        #else:
            plt.show()
    def find_abs(self,path):
        # assign directory
        #path = r'C:\Users\kdey\OneDrive - Werfen\Desktop\Python code for data analysis\test'
    
        abs = []
        aby = []
        # iterate over files in
        # that directory
        for filename in os.listdir(path):
            f = os.path.join(path, filename)
            
            df1=pd.read_csv(f, skiprows=2)
            df1.columns = ['Wavenumber','Absorbance']
            df1.head()
            df2 = df1.loc[(df1['Wavenumber'] > 1635) & (df1['Wavenumber'] < 1645)]
            
            df3=df2['Absorbance'].mean()
            #print(filename)
            #print(df2)
            #print(df3)
            abs.append(df3)
            aby.append(filename)
    
    
        print(aby)
        print(abs)
        x = pd.DataFrame({'filename':aby,'mean':abs})
        x.to_csv('find.csv',index=False, encoding='utf-8')

    def linear_reg(self,path):
        
        '''
        Importing data for fitting by linear regression
        
        filePATH = file path to CSV folder containing data for fitting r'filePATH'
        
        xlabel = x axis label
        
        ylabel = y axis label
        
        legend = 0,1, or 2 for legend to be located upper-right, upper-center, or upper-left, respectively
        
        saveyn = 'yes' to save plots, 'no' or anything else to not save plots - string
        
        save_csv_yn - 'yes' to save unit area normalized spectra as their own CSV, 'no' to not
        '''
        legend=0
        
        
        try:
            ## Import x,y data for fitting, z is error column
            x,y,z = np.loadtxt(path, delimiter=',', unpack=True, skiprows = 3)
        except:
            ### Import x,y data for fitting if no z (error column) is given
            x,y = np.loadtxt(path, delimiter=',', unpack=True, skiprows = 3)
        
        ### Fitting linear regression model to x,y data
        slope, intercept, r, p, std_err= stats.linregress(x,y)
        
        ### Linear function for fitting to data
        def linear(x):
            return slope * x + intercept
        
        ### Model object for plotting against x,y data
        mymodel = list(map(linear, x))
        
        ### Print fit parameters to double check
        print(slope, intercept, r, p, std_err)
        
        try:
            ### Try to plot x,y with y error (from z column) as scatter plot
            plt.errorbar(x,y,yerr = z, marker = 's',label = '1640', ls = 'none')
        except:
            ### If no y error given (from z column) just plot x,y as sactter with no error bars
            plt.scatter(x,y, label = '1640')
        
        ### Plot model as line together with scatter data
        plt.plot(x,mymodel,'r', label = 'Model')
        
        ### X-axis label
        plt.xlabel('Concentration(mM)', fontsize = 15)
        ### Y-axis label
        plt.ylabel('Absorbance - baseline subtracted', fontsize = 15)
        ### Change tick label (number) size
        plt.tick_params(labelsize = 12)
        ### Plot title
        plt.title('Linear Regression', fontsize = 15)
        ### Call tight_layout so all axis labels are visible in the saved image
        plt.tight_layout()
        
        ### Options for legend position
        #if legend == 0:
            #plt.legend(loc='upper right')
        #elif legend == 1:
            #plt.legend(loc='upper center')
        #elif legend == 2:
            #plt.legend(loc='upper left')
        
        ### Determine whether to just show plot or save it based on user input
        #if save_plot_yn == 'yes':
            #plt.tight_layout()
            #plt.savefig(filePATH.replace('.csv', '_linreg.jpeg'), dpi = 300)
        #else:
        plt.show()
        
        ### List of fit parameter names for tabulating fit results
        param_names = ['slope','intercept','r','p','std_err','intercept_stderr']
        
        ### Make empty list to fill with fit params
        param_list = []
        
        ### Append slope, intercept, and fit metris to param_list
        param_list.append(slope)
        param_list.append(intercept)
        param_list.append(r)
        param_list.append(p)
        param_list.append(std_err)
        
        ### Get folder path to save to from filePATH
        folder_path = os.path.dirname(os.path.abspath(path))
        
        #if save_csv_yn == 'yes':
            ### Write file names and absorbances at specified wavenumber to its own CSV file
            #with open(folder_path+'\\' + 'linear_reg_stats.CSV', 'w',newline = '') as new_file:
                #writer = csv.writer(new_file)
                #writer.writerows(zip(param_names, param_list)) 
        #else:
            #pass
        

if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    MainWindow = QtWidgets.QMainWindow()
    ui = Ui_MainWindow()
    ui.setupUi(MainWindow)
    MainWindow.show()
    sys.exit(app.exec_())

