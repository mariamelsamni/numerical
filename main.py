import sys

import numpy as np
from PyQt5.QtWidgets import QApplication
from PyQt5.QtWidgets import QPushButton
from PyQt5.QtWidgets import QVBoxLayout
from PyQt5.QtWidgets import QWidget,QHeaderView,QRadioButton
from PyQt5.QtWidgets import QLineEdit,QTableView
from PyQt5.QtWidgets import QGridLayout,QComboBox
from PyQt5.QtWidgets import QLabel,QStackedWidget
from PyQt5.QtCore import Qt,QAbstractTableModel
from PyQt5.QtGui import QDoubleValidator
from guass_jordan import Gauss_jordan as GJ
import pandas as pd
import Jacobi_seidel as Js
import GaussElimination as GE
import LU__doolittle
import Cholesky as CD
import Crout3 as Cr
import doolittle_coefficients
app = QApplication([])
window1=QWidget()
widget=QStackedWidget()
widget.addWidget(window1)
guess=[]
buffer=[]
time=None
class Model_Editable(QAbstractTableModel):
    def __init__(self, data):
        super().__init__()
        self._data = data

    def rowCount(self, index):
        return self._data.shape[0]

    def columnCount(self, index):
        return self._data.shape[1]


    def data(self, index, role=Qt.DisplayRole):
        if index.isValid():
            if role == Qt.DisplayRole or role == Qt.EditRole:
                value = self._data.iloc[index.row(), index.column()]
                return str(value)
    def setData(self, index, value, role):
        if role == Qt.EditRole:
            self._data.iloc[index.row()][index.column()] = value
            # print(self._data)
            return True
        return False

    def flags(self, index):
        return Qt.ItemIsSelectable | Qt.ItemIsEnabled | Qt.ItemIsEditable

    def headerData(self, section, orientation, role):
        if role == Qt.DisplayRole:
            if orientation == Qt.Horizontal:
                return str(self._data.columns[section])

            if orientation == Qt.Vertical:
                return str(self._data.index[section])

    # def setHeaderData(self, section, orientation, data, role=Qt.EditRole):
    #     if orientation == Qt.Horizontal and role in (Qt.DisplayRole, Qt.EditRole):
    #         try:
    #             self.horizontalHeaders[section] = data
    #             return True
    #         except:
    #             return False
    #     return super().setHeaderData(section, orientation, data, role)

    # def headerData(self, section, orientation, role=Qt.DisplayRole):
    #     if orientation == Qt.Horizontal and role == Qt.DisplayRole:
    #         try:
    #             return self.horizontalHeaders[section]
    #         except:
    #             pass
    #     return super().headerData(section, orientation, role)
class Model_nonEditable(QAbstractTableModel):
    def __init__(self, data,precision):
        super().__init__()
        self._data = data
        self._precision="%."+str(precision)+"f"
    def rowCount(self, index):
        return self._data.shape[0]

    def columnCount(self, index):
        return self._data.shape[1]


    def data(self, index, role):
        if index.isValid():
         if role == Qt.DisplayRole:
            value = self._data.iloc[index.row(), index.column()]
            if isinstance(value, float):

                return self._precision % value
            if isinstance(value, str):

                return '"%s"' % value
            return value


    def headerData(self, section, orientation, role):

        if role == Qt.DisplayRole:
            if orientation == Qt.Horizontal:
                return self._data.columns[section]

            if orientation == Qt.Vertical:
                return self._data.index[section]
def LU_window(Inputs, buffer,matrixsize, scale, precision,Scaling):
    buffer=Inputs.copy()
    widget.currentWidget().close()
    window4 = QWidget()
    options = QComboBox()
    characters_checker=[[0]*(matrixsize+1) for i in range(matrixsize)]
    try:
     for i in range(matrixsize):
        for j in range(matrixsize):
            characters_checker[i][j] = eval(Inputs.iloc[i][j])

     label = QLabel("Choose a LU Decomposition Method:")
     options.addItems(["Doolittle","Crout","Cholesky"])
    except:
        label = QLabel("Only Doolittle Method Supports Characters Solution:")
        options.addItems(["Doolittle"])

    GoToRes = QPushButton("Go To results")
    layout=QVBoxLayout()
    layout.addWidget(label)
    layout.addWidget(options)
    layout.addWidget(GoToRes)
    window4.setLayout(layout)
    widget.addWidget(window4)
    widget.setCurrentWidget(window4)
    widget.currentWidget().show()
    GoToRes.clicked.connect(lambda: ShowResultsWindow(Inputs,matrixsize,options.currentText(),guess,scale,precision,None,None,None,Scaling))
def Jacobi_Seidel_window(Inputs,buffer,matrixsize,name,precision):
    buffer = Inputs.copy()
    widget.currentWidget().close()
    rows = []
    rows.clear()
    for i in range(matrixsize):
        x = "X" + str(i + 1) + " ="
        rows.append(x)
    window4=QWidget()
    label=QLabel("Enter an Initial Guess :")
    layout=QVBoxLayout()
    fguess=["1"]*matrixsize

    fguess = pd.DataFrame(fguess, columns=['Guess'], index=rows)
    table = QTableView()
    model = Model_Editable(fguess)
    table.setModel(model)
    layout.addWidget(label)
    layout.addWidget(table)
    layout_options=QGridLayout()

    rearrange_options=QComboBox()
    rearrange_options.addItems(["Don't Re-arrange Matrix","Re-arrange Matrix"])
    label4=QLabel("Stopping conditions :")
    label2=QLabel("Enter Number of Iterations ")
    label3=QLabel("Enter The Absolute Relative Error")

    line_edit_Error=QLineEdit()
    line_edit_Itirations=QLineEdit()
    line_edit_Itirations.setPlaceholderText("Maximum 100")
    line_edit_Error.setPlaceholderText("0")
    layout_options.addWidget(label4,0,0)
    layout_options.addWidget(label2, 1, 0)
    layout_options.addWidget(line_edit_Itirations, 1, 1)
    layout_options.addWidget(label3, 2, 0)
    layout_options.addWidget(line_edit_Error, 2, 1)
    layout.addLayout(layout_options)
    GoToRes=QPushButton("Go To results")
    retrn = QPushButton("Return")
    layout.addWidget(rearrange_options)
    layout.addWidget(GoToRes)
    layout.addWidget(retrn)
    window4.setLayout(layout)
    widget.addWidget(window4)
    widget.setCurrentWidget(window4)
    widget.currentWidget().show()
    rule_Itirations = QDoubleValidator(0, 100,0)
    rule_Error = QDoubleValidator(0,5,10)
    line_edit_Itirations.setValidator(rule_Itirations)
    line_edit_Error.setValidator(rule_Error)
    GoToRes.clicked.connect(lambda :ShowResultsWindow(Inputs,matrixsize,name,fguess,False,precision,float(setDefaultError(line_edit_Error.text())),int(setDefaultItirations(line_edit_Itirations.text())),rearrange_options.currentText(),"Unscaled"))
    retrn.clicked.connect(lambda : getInputsWindow(matrixsize,buffer,True))
def setDefaultError(Error):
    try:
        x = eval(Error)
        if (x > 100):
            return 100
        return Error
    except:
        return -1
def setDefaultItirations(Itirations):
    try:
        x = eval(Itirations)
        if (x > 100):
            return 100
        return Itirations
    except:
        return 100
def moveTowindow():
    widget.setCurrentWidget(window1)


def ShowResultsWindow(Inputs,matrixsize,method_name,fguess,scale,precision,tolerance,numberOfItirations,Rearrange_condition,Scaling):
    characters=False
    buffer=Inputs.copy()
    coeifficients = [[0] * matrixsize for i in range(matrixsize)]
    equationsResults = [0] * matrixsize
    try:
     for i in range(matrixsize):
        for j in range(matrixsize):
            coeifficients[i][j] = eval(Inputs.iloc[i][j])
     for i in range(matrixsize):
        equationsResults[i] = eval(Inputs.iloc[i][matrixsize])
    except:
        characters=True
        for i in range(matrixsize):
            for j in range(matrixsize):
                coeifficients[i][j] = Inputs.iloc[i][j]
        for i in range(matrixsize):
            equationsResults[i] = Inputs.iloc[i][matrixsize]
    if(characters and (method_name!="Gauss Jordan" and method_name!="Doolittle")) :
         return
    rows = []
    rows.clear()
    for i in range(matrixsize):
        x = "X" + str(i + 1) + " ="
        rows.append(x)

    if (method_name == "Jacobi"):
         guess = [1] * matrixsize
         try:
          for i in range(matrixsize):
            guess[i] = eval(fguess.iloc[i][0])



          x,time= Js.Jacobi(coeifficients, equationsResults, guess, precision,Rearrange_condition,tolerance,numberOfItirations)

          flag = x.pop()

         except:
             return
    elif (method_name == "Seidel"):
         guess = [1] * matrixsize
         try:
          for i in range(matrixsize):
            guess[i] = eval(fguess.iloc[i][0])


          x,time= Js.Seidel(coeifficients, equationsResults, guess, precision, Rearrange_condition,tolerance,numberOfItirations)
          flag = x.pop()
         except:
             return

    elif (method_name == "Gauss Jordan"):
        if(not characters):
         x,time= GJ.guass_jordan_numbers(coeifficients, equationsResults,scale,precision)
         flag =x.pop()
        else:
            x,time=GJ.guass_jordan_characters(coeifficients,equationsResults)
            print("from gauss jordan x=",x)
            flag =x.pop()
    elif (method_name == "Gauss Elimination"):
         x,time=GE.GuassElimination(coeifficients,equationsResults,scale,precision)
         flag=x.pop()
    elif (method_name == "Doolittle"):
         if(not characters):
          x,flag,time=LU__doolittle.LU_doolittle(coeifficients,equationsResults,0.00005,matrixsize,precision,scale)
          print("elflag",flag)
          if flag=="infinite":
            flag="doolittle_infinite"
         else:
             x, flag,time= doolittle_coefficients.LU_doolittle(coeifficients, equationsResults, matrixsize)
             
    elif (method_name == "Crout"):
         x,time=Cr.crout(coeifficients,equationsResults,precision)
         flag = x.pop()
    elif (method_name == "Cholesky"):
         x,time=CD.Cholesky_Decomposition(coeifficients,matrixsize,equationsResults,precision)
         flag = x.pop()
    window3 = QWidget()
    layout = QVBoxLayout()
    rtrn = QPushButton("return")
    mainlabel = QLabel(method_name + " Method results :")
    layout.addWidget(mainlabel)

    if(flag=="correct" or flag=="infinite" or flag=="doolittle_infinite"):
     
     if (flag == "correct"):
      x=pd.DataFrame(x,columns=["Results"],index=rows)
     else:
         inf=QLabel("Infinite number of Solutions")
         layout.addWidget(inf)
         x = pd.DataFrame(x, columns=["One Of the Infinite Solutions"], index=rows)
     table = QTableView()
     model = Model_nonEditable(x,precision)
     table.setModel(model)
     header = table.horizontalHeader()
     header.setSectionResizeMode(QHeaderView.ResizeToContents)
     if (flag!="doolittle_infinite"):
        layout.addWidget(table)
     time="Runtime :"+str(np.round_(time,6))+" Seconds"
     Time_label=QLabel(time)
     layout.addWidget(Time_label)
    else:
        label1=QLabel(flag)
        layout.addWidget(label1)
    layout.addWidget(rtrn)
    window3.setLayout(layout)
    widget.addWidget(window3)
    widget.setCurrentWidget(window3)
    widget.show()
    rtrn.clicked.connect(lambda : decide(Inputs,buffer,matrixsize,method_name,Scaling,precision,True))
def getInputsWindow(matrixSize,buffer,edited):
    widget.setCurrentIndex(widget.currentIndex()+1)
    window2=QWidget()
    layout2 = QGridLayout()
    columns=[]
    columns.clear()
    rows=[]
    rows.clear()
    for i in range(matrixSize):
        eq="EQ "+str(i+1)
        rows.append(eq)
    for i in range(matrixSize+1):
      if(i!=matrixSize):
         x="X"+str(i+1)
      else:
          x="="
      columns.append(x)
    if(not edited):
         Inputs = pd.DataFrame([["0.0"] * (matrixSize + 1) for i in range(matrixSize)], columns=columns, index=rows)
    if(edited):
        Inputs = pd.DataFrame(buffer, columns=columns, index=rows)
    table=QTableView()
    model=Model_Editable(Inputs)
    table.setModel(model)
    layout2.addWidget(table)
    layout3=QGridLayout()
    Options=QComboBox()
    Options.addItems(["Gauss Elimination","Gauss Jordan","LU Decomposition","Jacobi","Seidel"])
    Scaling=QComboBox()
    Scaling.addItems(["Unscaled","Scaled"])
    res=QPushButton("Show Results")
    rtrn=QPushButton("Return")
    label2=QLabel("Choose a method :")
    label_scaling=QLabel("Choose Scaling Option :")
    label_precision=QLabel("Enter Precision :")
    label_max=QLabel(" maximum = 10")
    rule = QDoubleValidator(0, 10, 0)
    precision=QLineEdit("5")
    precision.setValidator(rule)
    layout3.addWidget(label2,0,0)
    layout3.addWidget(Options,0,1)
    layout3.addWidget(label_scaling,1,0)
    layout3.addWidget(Scaling,1,1)
    layout3.addWidget(label_precision,4,0)
    layout3.addWidget(precision,4,1)
    layout3.addWidget(label_max, 4, 2)
    layout3.addWidget(res,5,1)
    layout3.addWidget(rtrn,6,1)
    layout=QVBoxLayout()
    layout.addLayout(layout2)
    layout.addLayout(layout3)
    window2.setLayout(layout)
    widget.addWidget(window2)
    widget.setCurrentWidget(window2)
    widget.show()
    res.clicked.connect(lambda :decide(Inputs,buffer,matrixSize,Options.currentText(),Scaling.currentText(),int(float(SetDefault(precision.text()))),False))
    rtrn.clicked.connect(moveTowindow)
def decide(Inputs,buffer,matrixsize,name,Scaling,precision,Back):
    if(Scaling=="Scaled"):
        scale=True
    else:
        scale=False
    if((name=="Jacobi" or name=="Seidel")):
        Jacobi_Seidel_window(Inputs,buffer, matrixsize,name,precision)
    # elif ((name == "Jacobi" or name == "Seidel") and  Back):
    #     getInputsWindow(matrixsize, buffer, True)
    elif((name=="Gauss Elimination" or name=="Gauss Jordan") and not Back):
        ShowResultsWindow(Inputs, matrixsize, name,guess,scale,precision,None,None,None,Scaling)
    elif ((name == "Gauss Elimination" or name == "Gauss Jordan" or name=="Doolittle" or name=="Crout" or name=="Cholesky") and Back):
        getInputsWindow(matrixsize, buffer, True)
    elif(name=="LU Decomposition"):
        LU_window(Inputs,buffer, matrixsize, scale, precision,Scaling)
def SetDefault(text):
    try:
        x=eval(text)
        if(x>10):
            return 10
        return text
    except:
        return 5
def StartWindow() :
    layout1 = QGridLayout()
    label1 = QLabel("Enter number of Variables\nminmum = 2")
    label1.setAlignment(Qt.AlignCenter)
    rule = QDoubleValidator(2, 100,0)
    line_edit1 = QLineEdit()
    line_edit1.setPlaceholderText("Minimum 2")
    line_edit1.setValidator(rule)
    line_edit1.setAlignment(Qt.AlignCenter)
    layout1.addWidget(label1)
    layout1.addWidget(line_edit1)
    window1.setLayout(layout1)
    line_edit1.returnPressed.connect(lambda : getInputsWindow(int(float(line_edit1.text())),buffer,False))
    widget.show()
#Start Prpgram
StartWindow()
sys.exit(app.exec_())



