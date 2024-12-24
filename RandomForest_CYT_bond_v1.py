# Importing the libraries
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
##import sys
##sys.stdout = open('C:/Students/Spectra_bond_project/Logout_2_11_2023', 'w')
#import pandas as pd 
#import numpy as np
# データを読み込んでデータフレーム形式として保存
df = pd.read_csv('C:\Students\Spectra_bond_project\spectra_1_AuNP_CYT_bonds_MD_BL_corr.csv', delimiter="\t", index_col=0)
#df = pd.read_csv('C:\Students\Spectra_bond_project\s14_1_bc_GUA.csv', delimiter="\t", index_col=0)
#$df = pd.read_csv('C:\Students\Spectra_bond_project\s3_1_CYT_5MC_HMC_peak_freqency_calculation.csv', delimiter="\t", index_col=0)
print(df)
print("df printed:::::")

#読み込みデータをoutファイルに出力して確認
#df.to_csv('C:\Students\Spectra_bond_project\logout_RF_GUA_ini.out')   ##('file2.csv', header=False, index=False)
df.to_csv('C:\Students\Spectra_bond_project\logout_RF_CYT_1NP_ini.out')   ##('file2.csv', header=False, index=False)
##print(df)

#print('click ENTER to continue calculations after pause !!!!!')
#input()

f1=open('C:\Students\Spectra_bond_project\logout_RF_CYT_1NP_only.out', 'w+')
#f1=open('C:\Students\Spectra_bond_project\logout_RF_GUA_only.out', 'w+')
## with open('logout_RF_CYT__2.out', 'w+') as file2:
f3=open('C:\Students\Spectra_bond_project\out_RF_CYT_1NP_ran_ADE.out', 'w+')
f4=open('C:\Students\Spectra_bond_project\out_RF_CYT_1NP_fix_ADE.out', 'w+')
ferr=open('C:\Students\Spectra_bond_project\error_check.out', 'w+')
f1.write('----logout data------\n')

df.head()
# 説明変数
Xdf = df.iloc[:, 0:8].values
print("Xdf allocated*******")
print("Xdf ", Xdf[10:12,:])
X = df.to_numpy()
print("X allocated*******")
print("X ", X[10:12,:])

#X.shape
Xalp = df.iloc[:, 0:8].values
Xall = df.iloc[:, 0:8].values
#X = df.iloc[:, :2]
# 最大値を1、最小値を0にするような正規化。
# axis=1 とすれば、列ではなく行単位で正規化します。
lX=len(X)
lXm1=lX-1
print(len(X),'=length X')
L1=[str(len(X)),'=length X']
#X.shape
f1.writelines(L1)

#print('click ENTER to continue calculations after pause !!!!!')
#input()

alp_l=[0.0177, 0.0224, 0.0224, 0.0177, 0.0177, 0.0177]
alp_p=[0.0058, 0.0058, 0.0065, 0.0064, 0.0050, 0.0038]

lXb=len(alp_l)   #  number of bonds
lXbp1=lXb+1     #  average
print("alp_l length = ", lXb)

# 結果を格納するリスト
alp_t = []

# alp_tを計算してリストに追加（少数第4位までに丸める）
for l, p in zip(alp_l, alp_p):
    alp_t.append(round((l + 2 * p) / 3, 4))

# 結果を表示
print("alp_l=")
print(alp_l)
print("alp_p=")
print(alp_p)
print("alp_t=")
print(alp_t)

alp_mol = 0.0
for t in range (lXb):
    alp_mol = sum(alp_t)
    
print("alp_mol = ",alp_mol)
#alp_t [lXbp1]= alp_mol 
alp_t.append(alp_mol)   
print("appended alp_t=")
print(alp_t)

print("X  *  alpha  ======")

###Xalp = np.random.rand(12, 3)  # Sample data, 12 rows and 3 columns
Xav = 0
ri = 0
for r in  X:      #X[0:12,:] :    #(1,lX):  #(1,12): #test
    ri += ri
    if ri <=15:
        print(r)
#    for value in r:
#        ferr.write(r"  {value:.4f}\n")
#f3.write("end__x\n")

    Xav = 0
    for c in range (1,lXbp1):
        cm1 = c-1
        cp1 = c+1
        print(ri,c,cm1)
#        print(X[r][c])
#        print(al[p_t[c])
        Xrc = r[c]
        Xav = Xav + Xrc*alp_t[cm1]   # spectral average with polarizability weights
###        print(Xav,Xrc,alp_t[cm1])
        Xalp[ri][c] = Xrc*alp_t[cm1]   # spectral average with polarizability weights
###        print(Xav,Xrc,alp_t[cm1], Xalp[ri][c])
    Xalp[ri][0] = r[0]
    Xalp[ri][lXbp1] = Xav
    r[lXbp1]= Xav   # spectral average with polarizability weights
    #print(r[lXbp1],Xalp[ri][lXbp1])
    if ri <=15:
        print(r)
        print(Xalp[ri][:])
#    print(r)
#    print(Xalp[ri][:])

#    for value in r:
#        ferr.write(r"  {value:.4f}\n")
#        #ferr.write(Xalp[ri][:]"  {value:.4f}\n")

X_ini = X   # save updated with last column X 
print("Xalp  ======")
ri = 0
for r1 in  Xalp:      #X[0:12,:] :    #(1,lX):  #(1,12): #test
    ri += ri
    if ri <= 15:
        print(r1)
X = Xalp    # save all bond spectra with Alpha_bond coefficients into X
print("X = Xalp  ======")
ri = 0
for r2 in  X:      #X[0:12,:] :    #(1,lX):  #(1,12): #test
    ri += ri
    if ri <= 15:
        print(r2)

##X = [row[:] for row in Xalp]
#X_ini = np.array(X)   # save updated with last column X 
#X = np.array(Xalp)    # save all bond spectra with Alpha_bond coefficients into X

#------------------------------- make input matrix for fixed split ------
print("X  split  =====")
Xfix = []
Xdf1 = []
rj = 0
rf = 0
rn = 0
for r in X:    #(1,lX):  #(1,12): #test
    rj += 1
    if rj % 5 == 0:  # Condition for 20% test data removal
        rf += 1
        Xfix.append(r)  # Add row to Xfix (for rows where rj % 5 == 0)
    else:
        if rj <= 15:
            print(r)
        rn += 1
        Xdf1.append(r)  # Add row to Xdf1 (for other rows)


#    if rj == -1:
#       Xdf1.append(r) 
    rj += 1
    if rj  % 5 == 0:      #  20% test data removal from the whole data file
        rf += 1
        Xfix.append(r)
    else:
        if rj == 1:
            print(r)
            Xdf1.append(r) 

        if rj > 15:
            print(r)

            rn += 1
            Xdf1.append(r)
Xfix = np.array(Xfix)
Xdf1 = np.array(Xdf1)

Xf = np.concatenate((Xdf1, Xfix)) #  made X matrix with 20% fixed test data at the end of file 

lXfix=len(Xfix)       
lXdf1=len(Xdf1)       
print("---lengthes of Xfix & Xdf1 ", lXfix, lXdf1)
print(Xfix)
print("Xdf1==")
print(Xdf1)

#print('click ENTER to continue calculations after pause !!!!!')
#input()    # pause

##print ('===averaged spectra X')
##print (X[10:20,lXbp1])
##print ('===averaged spectra Xalp')
##print (Xalp[10:20,0:lXbp1])

print('click ENTER to continue calculations after pause !!!!!')
input()    # pause

#-------------------------------
#Xf = df.iloc[:, j:(j+1)].values for j in range (7)  # error
X1 = df.iloc[:, 0:1].values
X2 = df.iloc[:, 0:1].values
X3 = df.iloc[:, 0:1].values
X4 = df.iloc[:, 0:1].values
X5 = df.iloc[:, 0:1].values
X6 = df.iloc[:, 0:1].values
X7 = df.iloc[:, 0:1].values
for i in range(lXm1):
    X1[i,0:1] = 0
    X2[i,0:1] = 0
    X3[i,0:1] = 0
    X4[i,0:1] = 0
    X5[i,0:1] = 0
    X6[i,0:1] = 0
    X7[i,0:1] = 0
#---------------------------------
    
# 目的変数
y = df.iloc[:, -8].values
#y = df.iloc[:, :2]
###data = pd.read_csv('C:\Students\Spectra_bond_project\Salaries.csv')

##X_str=X.tostring()
##y_str=y.tostring()  ### gives integer for writelines
X_str=np.array2string(X,formatter={'float_kind':lambda X: "%.4f" % X})
X_str=np.array2string(X,formatter={'float_kind':lambda X: "%.4f" % X})
y_str=np.array2string(y,formatter={'float_kind':lambda y: "%.4f" % y})

print("---X")
print(X)
#print("---Xalp")
#print(Xalp)

print('click ENTER to continue calculations after pause !!!!!')
input()    # pause


f1.write("\n---x\n")
f1.writelines(X_str)    
##print(x[:, 0])
##print(x[:, 1])
print("---y")
print(y)
f1.write("\n---y\n")
f1.writelines(y_str)    

# 訓練データとテストデータに分割するメソッドのインポート
from sklearn.model_selection import train_test_split 
# 訓練データ・テストデータへ6:4の比でランダムに分割
##X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2) 
##random_state: Pass an int for reproducible output across multiple function call
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42) 

X_train_str=np.array2string(X_train,formatter={'float_kind':lambda X_train: "%.4f" % X_train})
y_train_str=np.array2string(y_train,formatter={'float_kind':lambda y_train: "%.4f" % y_train})
X_test_str=np.array2string(X_test,formatter={'float_kind':lambda X_test: "%.4f" %  X_test})
y_test_str=np.array2string(y_test,formatter={'float_kind':lambda y_test: "%.4f" % y_test})

print("------x_train")
print(X_train)
print("------x_test")
print(X_test)
print("------y_train")
print(y_train)
print("------y_test")
print(y_test)

print('click ENTER to continue calculations after pause !!!!!')
input()    # pause


f1.write("\n\n")
lX_test=len(X_test)
lXm1_test=lX_test-1
print(len(X_test),'=length X_test')
L1_test=[str(len(X_test)),'=length X_test']
f1.writelines(L1_test)


f1.write("\n---X_train\n")
f1.writelines(X_train_str)  
f1.write("\n---X_test\n")
f1.writelines(X_test_str)  
f1.write("\n---y_train\n")
f1.writelines(y_train_str)  
f1.write("\n---y_test\n")
f1.writelines(y_test_str)  
f1.write("\n\n")
#### Random Forest 
print("Start Random Forest------random_test")
# Fitting Random Forest Regression to the dataset
# import the regressor
from sklearn.ensemble import RandomForestRegressor
  
 # create regressor object
regressor = RandomForestRegressor(n_estimators = 5, random_state = 0)
  
# fit the regressor with x and y data
##regressor.fit(X, y) 
##Y_pred_RF = regressor.predict(np.array([6.5]).reshape(1, 1))  # test the output by changing values

#regressor = RandomForestRegressor(n_estimators=20, random_state=0)
regressor.fit(X_train, y_train)
Y_pred_RF = regressor.predict(X_test)  # test the output by changing values


print("Random Forest Error")

from sklearn import metrics
print('Mean Absolute Error:', metrics.mean_absolute_error(y_test, Y_pred_RF))
print('Mean Squared Error:', metrics.mean_squared_error(y_test, Y_pred_RF))
print('Root Mean Squared Error:', np.sqrt(metrics.mean_squared_error(y_test, Y_pred_RF)))

f1.write("Random Forest Error\n")
L_MAE=['Mean Absolute Error:',str(metrics.mean_absolute_error(y_test, Y_pred_RF)),'\n']
L_MSE=['Mean Squared Error:',str(metrics.mean_squared_error(y_test, Y_pred_RF)),'\n']
L_RMSE=['Root Mean Squared Error:',str(np.sqrt(metrics.mean_squared_error(y_test, Y_pred_RF))),'\n']

f1.writelines(L_MAE)
f1.writelines(L_MSE)    
f1.writelines(L_RMSE) 

print('click ENTER to continue calculations after pause !!!!!')
input()    # pause
   

# Visualising the Random Forest Regression results
  
# arrange for creating a range of values
# from min value of x to max 
# value of x with a difference of 0.01 
# between two consecutive values
##X_grid = np.arange(min(X), max(X), 0.01) 
##print(np.arange(3, 10, 2))
#x_min=min(X_train)
#print(type(x_min),x_min)
#x_max=max(X_train)
#print(type(x_max),x_max)
##x_min=all(min(X[[:,0]]))
###-----x_step=int(X[1,0]-X[0,0])  # not correct
###-----print(x_step,"=x_step")

x_min=min(X[:,0])   ##x_min=min(X[:,1])
print(type(x_min),x_min)
x_max=max(X[:,0])   ##x_max=max(X[:,1])
print(type(x_max),x_max)

x_step=int((x_max-x_min)/lXm1)    ### X_grid[0:85]
print(x_step,"=x_step_len")

#X_grid = np.arange(x_min, x_max, 0.1)   # for single vector
X_grid = np.arange(x_min, x_max, x_step)   # for single vector
print(len(X_grid),"=len(X_grid)")
print(X_grid)
  
# reshape for reshaping the data into a len(X_grid)*1 array, 
# i.e. to make a column out of the X_grid value                  
#X_grid = X_grid.reshape((len(X_grid), 1))
X_grid = X_grid.reshape((len(X_grid), 1))
  
# Scatter plot for original data
##plt.scatter(X, y, color = 'blue')
##plt.scatter(X[:,0], y, color = 'blue')    # correct
peak_type=Xall[:,4]
print('peak_type')
print(peak_type)
j=0
for i in range(lXm1):
##if peak_type.all() == 1.0000e+00:
##if X[:,0].any() == 1.0000e+00:
##     plt.scatter(X[:,0], y, color = 'blue')  
    if Xall[i,4] == 1.0000e+00:
        X1[i,0:1]=X[i,0:1]
##        X1[j,1]=X[i,1]
        print(X1[i,0:1],X[i,0:1])
    if   X1[i,0:1] > [0]:
# plt.scatter(X1[:,0], y, color = 'blue')
        plt.scatter(X1[i,0], y[i], s=10, color = 'blue')  
print(X1,"=x1")

#print('click ENTER to continue calculations after pause !!!!!')
#input()

#if X1[i,0:1] > 0:   
#   plt.scatter(X1[:,0], y, color = 'blue')  
j=0 
for i in range(lXm1):
##if peak_type.all() == 2.0000e+00:
##    if X[i,0].any() == 2.0000e+00:
    if Xall[i,4] == 2.0000e+00:
        X2[i,0:1]=X[i,0:1]
    if  X2[i,0:1] > 0:
        plt.scatter(X2[i,0], y[i], s=10, color = 'yellow')
#       plt.scatter(X2[:,0], y, color = 'yellow')  
print(X2,"=x2")
#if not X2[i,0:1] == 0:
#plt.scatter(X2[:,0], y, color = 'yellow')  
j=0
for i in range(lXm1):
##if peak_type.all() == 3.0000e+00:
##    if X[i,0].any() == 3.0000e+00:
    if Xall[i,4] == 3.0000e+00:
        X3[i,0:1]=X[i,0:1]
    if  X3[i,0:1] > 0:
        plt.scatter(X3[i,0], y[i], s=10, color = 'red')
#       plt.scatter(X3[:,0], y, color = 'red')
print(X3,"=x3")
#if not X2[i,0:1] == 0:
#plt.scatter(X3[:,0], y, color = 'red')  
  
j=0
for i in range(lXm1):
##if peak_type.all() == 3.0000e+00:
##    if X[i,0].any() == 3.0000e+00:
    if Xall[i,4] >= 4.0000e+00:
        X4[i,0:1]=X[i,0:1]
    if  X4[i,0:1] > 0:
        plt.scatter(X4[i,0], y[i], s=10, color = 'black')
#       plt.scatter(X3[:,0], y, color = 'red')
print(X4,"=x4")
#if not X2[i,0:1] == 0:
#plt.scatter(X3[:,0], y, color = 'red')  

# plot predicted data
plt.plot(X_grid[0:lX], regressor.predict(X), 
         color = 'green')
##plt.plot(X_grid[0:lX], regressor.predict(X_test), 
##         color = 'blu')   ##  ValueError: x and y must have same first dimension, but have shapes (737, 1) and (148,)

##----------output to file---------------------
X_grid_vec=X_grid[0:lX]
X_grid_str=np.array2string(X_grid_vec,formatter={'float_kind':lambda X_grid_vec: "%.4f" % X_grid_vec})
Y_pred_RF_vec=regressor.predict(X)  ##Y_pred_RF   ##
Y_pred_RF_vec.shape
Y_pred_RF_vec_line=Y_pred_RF_vec[:,None]  ## numpy.array()  list to array    #Y_pred_RF_vec[0:lX]
Y_grid_str=np.array2string(Y_pred_RF_vec,formatter={'float_kind':lambda Y_pred_RF_vec: "%.4f" % Y_pred_RF_vec})
Y_grid_str_line=np.array2string(Y_pred_RF_vec_line,formatter={'float_kind':lambda Y_pred_RF_vec_line: "%.4f" % Y_pred_RF_vec_line})

#--- Y1_pred_RF_vec ---- length ------
f3.write("\n\n")
lY_test=len(Y_pred_RF_vec)
lym_test1=lY_test-1
print(len(Y_pred_RF_vec),'=length Y_pred_RF_vect')
L1_testY=[str(len(Y_pred_RF_vec)),'=length Y_pred_RF_vec']
f3.writelines(L1_testY)
f3.write("\n\n")

f3.write("Random Forest output data__X__random train\n")
f3.writelines(X_grid_str)
f3.write("end__x\n")
f3.write("Random Forest output data__y\n")
f3.writelines(Y_grid_str)
f3.write("\n\n")
f3.write("Random Forest output data in line__y\n")
f3.writelines(Y_grid_str_line)

#------------ !!! manual input of train_test split !!! --------------**------------**--------
df1 = pd.read_csv('C:\Students\Spectra_bond_project\s15_2_bc_ADE.csv', delimiter="\t", index_col=0)
#df1 = pd.read_csv('C:\Students\Spectra_bond_project\s14_2_bc_GUA.csv', delimiter="\t", index_col=0)
#df1 = pd.read_csv('C:\Students\Spectra_bond_project\s13_2_bc_ADE.csv', delimiter="\t", index_col=0)
#df1 = pd.read_csv('C:\Students\Spectra_bond_project\s12_2_bc_THY.csv', delimiter="\t", index_col=0)
#df1 = pd.read_csv('C:\Students\Spectra_bond_project\s11_2_bc_HMC.csv', delimiter="\t", index_col=0)
#df1 = pd.read_csv('C:\Students\Spectra_bond_project\s10_2_bc_5MC.csv', delimiter="\t", index_col=0)
#df1 = pd.read_csv('C:\Students\Spectra_bond_project\s9_2_bc_CYT.csv', delimiter="\t", index_col=0)df1 = pd.read_csv('C:\Students\Spectra_bond_project\s9_2_bc_CYT.csv', delimiter="\t", index_col=0)
#df1 = pd.read_csv('C:\Students\Spectra_bond_project\s8_2_bc_CYT_5MC_HMC_peak_frequency_calculations.csv', delimiter="\t", index_col=0)
#df1 = pd.read_csv('C:\Students\Spectra_bond_project\s7_2_bc_CYT_5MC_HMC_peak_frequency_calculations.csv', delimiter="\t", index_col=0)
#df1 = pd.read_csv('C:\Students\Spectra_bond_project\s5_2_CYT_5MC_HMC_peak_freqency_calculation.csv', delimiter="\t", index_col=0)
#df1 = pd.read_csv('C:\Students\Spectra_bond_project\s3_2_CYT_5MC_HMC_peak_freqency_calculation.csv', delimiter="\t", index_col=0)
print(df)
#読み込みデータをoutファイルに出力して確認
#df1.to_csv('C:\Students\Spectra_bond_project\logout_RF_CYT_ini_df1.out')   ##('file2.csv', header=False, index=False)
#df1.to_csv('C:\Students\Spectra_bond_project\logout_RF_THY_ini_df1.out')   ##('file2.csv', header=False, index=False)
#df1.to_csv('C:\Students\Spectra_bond_project\logout_RF_ADE_ini_df1.out')   ##('file2.csv', header=False, index=False)
#df1.to_csv('C:\Students\Spectra_bond_project\logout_RF_GUA_ini_df1.out')   ##('file2.csv', header=False, index=False)
df1.to_csv('C:\Students\Spectra_bond_project\logout_RF_ADE2_ini_df1.out')   ##('file2.csv', header=False, index=False)
df1.head()
# 説明変数
Xdf1 = df1.iloc[:, 0:4].values
Xalldf1 = df1.iloc[:, 0:5].values

lX1=len(Xdf1)
lX1m1=lX1-1
print(len(Xdf1),'=length Xdf1')
L1_1=[str(len(Xdf1)),'=length Xdf1']
f1.writelines(L1_1)

ydf1 = df1.iloc[:, -4].values

#print('click ENTER to continue calculations after pause !!!!!')
#input()

#-------------------------------
X1_1 = df1.iloc[:, 0:1].values
X2_1 = df1.iloc[:, 0:1].values
X3_1 = df1.iloc[:, 0:1].values
X4_1 = df1.iloc[:, 0:1].values
for i in range(lX1m1):
    X1_1[i,0:1] = 0
    X2_1[i,0:1] = 0
    X3_1[i,0:1] = 0
    X4_1[i,0:1] = 0
#---------------------------------    
print("---x")
print(Xdf1)
print("---x_all")
print(Xalldf1)

# 訓練データ・テストデータへ6:4の比でランダムに分割
#X1_train, X1_test, y1_train, y1_test = train_test_split(Xdf1, ydf1, test_size=0.2) 

#regressor.fit(X1_train, y1_train)
#Y1_pred_RF = regressor.predict(X1_test)  # test the output by changing values

X1_train, X1_test, y1_train, y1_test = train_test_split(Xdf1, ydf1, train_size=263, test_size=66, shuffle=False)
#X1_train, X1_test, y1_train, y1_test = train_test_split(Xdf1, ydf1, train_size=296, test_size=33, shuffle=False) 
#X1_train, X1_test, y1_train, y1_test = train_test_split(Xdf1, ydf1, train_size=221, test_size=24, shuffle=False) 
#X1_train, X1_test, y1_train, y1_test = train_test_split(Xdf1, ydf1, train_size=590, test_size=147, shuffle=False) 
#s4##X1_train, X1_test, y1_train, y1_test = train_test_split(Xdf1, ydf1, train_size=105, test_size=29, shuffle=False) 
#s4##X1_train, X1_test, y1_train, y1_test = train_test_split(Xdf1, ydf1, train_size=86, test_size=24, shuffle=False) 
#s3##X1_train, X1_test, y1_train, y1_test = train_test_split(Xdf1, ydf1, train_size=67, test_size=18, shuffle=False) 
print("------x1_train2")
print(X1_train)
print("------x1_test2")
print(X1_test)
print("------y1_train2")
print(y1_train)
print("------y1_test2")
print(y1_test)

X1_train_str=np.array2string(X1_train,formatter={'float_kind':lambda X1_train: "%.4f" % X1_train})
y1_train_str=np.array2string(y1_train,formatter={'float_kind':lambda y1_train: "%.4f" % y1_train})
X1_test_str=np.array2string(X1_test,formatter={'float_kind':lambda X1_test: "%.4f" %  X1_test})
y1_test_str=np.array2string(y1_test,formatter={'float_kind':lambda y1_test: "%.4f" % y1_test})

f1.write("\n\n")
lX_test1=len(X1_test)
lXm1_test1=lX_test1-1
print(len(X1_test),'=length X1_test')
L1_test1=[str(len(X1_test)),'=length X1_test']
f1.writelines(L1_test1)


f1.write("\n---X1_train_fix\n")
f1.writelines(X1_train_str)  
f1.write("\n---X1_test_fix\n")
f1.writelines(X1_test_str)  
f1.write("\n---y1_train_fix\n")
f1.writelines(y1_train_str)  
f1.write("\n---y1_test_fix\n")
f1.writelines(y1_test_str)  
f1.write("\n\n")

print("Start Random Forest------fixed_test")
#
regressor = RandomForestRegressor(n_estimators = 20, random_state = 2, max_leaf_nodes = 10)
regressor.fit(X1_train, y1_train)
Y1_pred_RF = regressor.predict(X1_test)  # test the output by changing values

print("Random Forest Error_no_shuffle")

#from sklearn import metrics
print('Mean Absolute Error:', metrics.mean_absolute_error(y1_test, Y1_pred_RF))
print('Mean Squared Error:', metrics.mean_squared_error(y1_test, Y1_pred_RF))
print('Root Mean Squared Error:', np.sqrt(metrics.mean_squared_error(y1_test, Y1_pred_RF)))

f1.write("Random Forest Error\n")
L1_MAE=['Mean Absolute Error:',str(metrics.mean_absolute_error(y1_test, Y1_pred_RF)),'\n']
L1_MSE=['Mean Squared Error:',str(metrics.mean_squared_error(y1_test, Y1_pred_RF)),'\n']
L1_RMSE=['Root Mean Squared Error:',str(np.sqrt(metrics.mean_squared_error(y1_test, Y1_pred_RF))),'\n']

f1.writelines(L1_MAE)
f1.writelines(L1_MSE)    
f1.writelines(L1_RMSE)    

x1_min=min(Xdf1[:,0])   ##x_min=min(X[:,1])
print(type(x1_min),x1_min)
x1_max=max(Xdf1[:,0])   ##x_max=max(X[:,1])
print(type(x1_max),x1_max)

x1_step=int((x1_max-x1_min)/lX1m1)    ### X_grid[0:85]
print(x1_step,"=x1_step_len")

X1_grid = np.arange(x1_min, x1_max, x1_step)   # for single vector
print(len(X1_grid),"=len(X1_grid)")
print(X1_grid)
  
# reshape for reshaping the data into a len(X_grid)*1 array, 
# i.e. to make a column out of the X_grid value                  
X1_grid = X1_grid.reshape((len(X1_grid), 1))



plt.plot(X1_grid[0:lX1], regressor.predict(Xdf1), 
         color = 'purple')

###Y_pred_RF_y = regressor.predict(y_test)  # test the output by changing values
###plt.plot(X_grid[0:85], regressor.predict(y), 
###         color = 'orangered')


plt.title('s15__1 vs 2__ADE_peak_frequency_calculation') 
#plt.title('s14__1 vs 2__GUA_peak_frequency_calculation')
#plt.title('s13__1 vs 2__ADE_peak_frequency_calculation')
#plt.title('s12__1 vs 2__THY_peak_frequency_calculation')
#plt.title('s11__1 vs 2__HMC_peak_frequency_calculation')
#plt.title('s9__1 vs 2__CYT_peak_frequency_calculation')
#plt.title('s5__1 vs 2__CYT_5MC_HMC_peak_frequency_calculation')
#plt.title('s3__1 vs 2__CYT_5MC_HMC_peak_frequency_calculation')
#plt.title('s4_CYT_5MC_HMC_peak_frequency_calculation')
plt.xlabel('frequency')
plt.ylabel('Peak Height')
plt.show()

##----------output to file---------------------
X1_grid_vec=X1_grid[0:lX]
X1_grid_str=np.array2string(X1_grid_vec,formatter={'float_kind':lambda X1_grid_vec: "%.4f" % X1_grid_vec})
Y1_pred_RF_vec=regressor.predict(Xdf1)  ##Y1_pred_RF   ##regressor.predict(X)
Y1_grid_str=np.array2string(Y1_pred_RF_vec,formatter={'float_kind':lambda Y1_pred_RF_vec: "%.4f" % Y1_pred_RF_vec})
Y1_pred_RF_vec_line=Y1_pred_RF_vec[:,None]  ## numpy.array()  list to array    #Y_pred_RF_vec[0:lX]
Y1_grid_str_line=np.array2string(Y1_pred_RF_vec_line,formatter={'float_kind':lambda Y1_pred_RF_vec_line: "%.4f" % Y1_pred_RF_vec_line})

#--- Y1_pred_RF_vec ---- length ------
f4.write("\n\n")
lY_test1=len(Y1_pred_RF_vec)
lym1_test1=lY_test1-1
print(len(Y1_pred_RF_vec),'=length Y1_pred_RF_vect')
L1_test1Y=[str(len(Y1_pred_RF_vec)),'=length Y1_pred_RF_vec']
f4.writelines(L1_test1Y)
f4.write("\n\n")

f4.write("Random Forest output data__X__fixed train\n")
f4.writelines(X1_grid_str)
f4.write("end__x\n")
f4.write("Random Forest output data__y\n")
f4.writelines(Y1_grid_str)
f4.write("\n\n")
f4.write("Random Forest output data in line__y\n")
f4.writelines(Y1_grid_str_line)


#--------------- Logistic regression -----------------------------------
#  ---This solver needs samples of at least 2 classes in the data----

XL = df.iloc[:, :4]
#XL = df.iloc[:, :2]
yl = df.iloc[:, -1]

print("---x_LR")
print(XL)
##print(x[:, 0])
##print(x[:, 1])
print("---y_LR")
print(yl)

# 訓練データとテストデータに分割するメソッドのインポート
from sklearn.model_selection import train_test_split 
# 訓練データ・テストデータへ6:4の比でランダムに分割
##X_LR_train, X_LR_test, y_LR_train, y_LR_test = train_test_split(XL, yl, test_size=0.4)
test_size_LR=0.25 
X_LR_train, X_LR_test, y_LR_train, y_LR_test = train_test_split(XL, yl, test_size=test_size_LR) 
#X_LR_train, X_LR_test, y_LR_train, y_LR_test = train_test_split(XL, yl, test_size=0.5) 



from sklearn.linear_model import LogisticRegression # ロジスティック回帰
classifier = LogisticRegression(max_iter=200) # 分類器の生成
classifier.fit(X_LR_train, y_LR_train) #学習
import timeit # 実行時間を計測するためのライブラリ
timeit.timeit(lambda: classifier.fit(X_LR_train, y_LR_train), number=1)
# 正解率 (train) : 学習に用いたデータをどのくらい正しく予測できるか
classifier.score(X_LR_train, y_LR_train)
# 正解率 (test) : 学習に用いなかったデータをどのくらい正しく予測できるか
classifier.score(X_LR_test,y_LR_test)
# 学習に用いなかったデータを予測する
y_pred = classifier.predict(X_LR_test)
print(y_pred)
from sklearn.metrics import confusion_matrix # 混同行列を計算するメソッド
# 予測結果と、正解（本当の答え）がどのくらい合っていたかを表す混同行列
pd.DataFrame(confusion_matrix(y_pred, y_LR_test), 
             index=['predicted 0', 'predicted 1','predicted 2'], columns=['real 0', 'real 1', 'real 2'])
# 予測結果の自信の強さを計算する
y_proba = classifier.predict_proba(X_LR_test)

print(y_proba)

print("LogisticRegression End")
# Error  y_true and y_pred have different number of output (1!=3)
#print('Mean Absolute Error:', metrics.mean_absolute_error(y_LR_test, y_proba))
#print('Mean Squared Error:', metrics.mean_squared_error(y_LR_test, y_proba))
#print('Root Mean Squared Error:', np.sqrt(metrics.mean_squared_error(y_LR_test, y_proba)))


#%matplotlib inline
# ROC curve を描く
from sklearn.metrics import roc_curve
from sklearn.metrics import auc

# AUCスコアを出す
#y_test_vec=y_LR_test
y_test_vec=pd.DataFrame(y_LR_test).to_numpy()
y_test=np.asarray(y_test_vec).reshape(-1)
###array([[1. , 3. ],
###       [2. , 4.5]])
y_proba_vec=y_proba[:, 2]
print("y_LR_test")
print(type(y_LR_test),y_LR_test)
print(len(y_LR_test),"=len(y_LR_test)")
print("y_test_vec")
print(type(y_test_vec),y_test_vec)
print(len(y_test_vec),"=len(y_test_vec)")
print("y_test")
print(type(y_test),y_test)
print(len(y_test),"=len(y_test)")

print("y_proba_vec")
print(len(y_proba_vec),"=len(y_proba_vec)")
print(type(y_proba_vec),y_proba_vec)
print("y_proba")
print(len(y_proba),"=len(y_proba)")
print(type(y_proba),y_proba)

#errors
print('Mean Absolute Error:', metrics.mean_absolute_error(y_LR_test,y_proba_vec))
print('Mean Squared Error:', metrics.mean_squared_error(y_LR_test,y_proba_vec))
print('Root Mean Squared Error:', np.sqrt(metrics.mean_squared_error(y_LR_test,y_proba_vec)))

y_in = np.array(y_test,)         # y_test values in the range [1-3] require pos_label=3 in 
y_scores = np.array(y_proba_vec)          #   roc_curve(y_in, y_scores,pos_label=3)
print(y_in)
print(y_scores)
fpr, tpr, thresholds = roc_curve(y_in, y_scores,pos_label=3)
##fpr, tpr, thresholds = roc_curve(y_test, y_proba_vec)
#fpr, tpr, thresholds = roc_curve(y_LR_test, y_proba_vec)
print(fpr)
print(tpr)   # nan
print(thresholds)
roc_auc = auc(fpr, tpr)
print ("AUC curve : %f" % roc_auc)

import matplotlib.pyplot as plt
#%matplotlib inline
# ROC curve を描く
plt.figure(figsize=(4,4))
plt.plot(fpr, tpr, label='ROC curve (area = %0.2f)' % roc_auc)
plt.plot(fpr, tpr, label='test_size =%0.2f' % test_size_LR)
plt.plot([0, 1], [0, 1], 'k--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.0])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('ROC curve: AUC=%0.2f' % roc_auc)
##plt.title('test_size =%0.2f' % test_size_LR)
plt.legend(loc="lower right")
plt.show()

f1.close()

##sys.stdout.close()