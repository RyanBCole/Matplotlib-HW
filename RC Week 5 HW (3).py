#!/usr/bin/env python
# coding: utf-8

# In[2]:


import matplotlib.pyplot as plt
import pandas as pd
import scipy.stats as st


# In[3]:


mouse_metadata = pd.read_csv(r'C:\Users\redye\Documents\Mouse_metadata.csv')
study_results = pd.read_csv(r'C:\Users\redye\Documents\Study_results.csv')


# In[4]:


mouse_metadata


# In[5]:


study_results


# In[6]:


mouse_comb = pd.merge(mouse_metadata, study_results, on="Mouse ID")
mouse_comb


# In[8]:


# Checking the number of mice.
mouse_comb['Mouse ID'].nunique()


# In[10]:


# Getting the duplicate mice by ID number that shows up for Mouse ID and Timepoint. 
dup_mice = mouse_comb.loc[mouse_comb.duplicated(subset=['Mouse ID', 'Timepoint']),'Mouse ID'].unique()
dup_mice


# In[12]:


# Optional: Get all the data for the duplicate mouse ID. 
dup_mice_df = mouse_comb.loc[mouse_comb["Mouse ID"] == "g989", :]
dup_mice_df


# In[14]:


# Create a clean DataFrame by dropping the duplicate mouse by its ID.
clean_df = mouse_comb[mouse_comb['Mouse ID'].isin(dup_mice)==False]
clean_df.head()


# In[15]:


# Checking the number of mice in the clean DataFrame.
mice = clean_df["Mouse ID"].nunique()

mice


# In[16]:


# Generate a summary statistics table of mean, median, variance, standard deviation, and SEM of the tumor volume for each regimen

# Use groupby and summary statistical methods to calculate the following properties of each drug regimen: 
# mean, median, variance, standard deviation, and SEM of the tumor volume. 
# Assemble the resulting series into a single summary DataFrame.

mean = clean_df['Tumor Volume (mm3)'].groupby(clean_df['Drug Regimen']).mean()
median = clean_df['Tumor Volume (mm3)'].groupby(clean_df['Drug Regimen']).median()
var = clean_df['Tumor Volume (mm3)'].groupby(clean_df['Drug Regimen']).var()
std = clean_df['Tumor Volume (mm3)'].groupby(clean_df['Drug Regimen']).std()
sem = clean_df['Tumor Volume (mm3)'].groupby(clean_df['Drug Regimen']).sem()

summary_df = pd.DataFrame({"Mean Tumor Volume":mean, 
                            "Median Tumor Volume":median, 
                           "Tumor Volume Variance":var, 
                           "Tumor Volume Std. Dev.":std, 
                           "Tumor Volume Std. Err.":sem})
summary_df


# In[17]:


# Generate a summary statistics table of mean, median, variance, standard deviation, 
# and SEM of the tumor volume for each regimen

# Using the aggregation method, produce the same summary statistics in a single line.

summary_total =  clean_df.groupby(['Drug Regimen'])[['Tumor Volume (mm3)']].agg(['mean', 'median', 'var', 'std', 'sem'])
summary_total


# In[29]:


drug_mice = clean_df["Drug Regimen"].value_counts()
drug_mice


# In[30]:


# Generate a bar plot showing the total number of timepoints for all mice tested for each drug regimen using Pandas.

plot_pandas = drug_mice.plot.bar(color='b')  

plt.xlabel("Drug Regimen")
plt.ylabel("Number of Mice")
plt.title("Number of Mice per Treatment")


# In[33]:


# Generate a bar plot showing the total number of timepoints for all mice tested for each drug regimen using matplotlib.

x_axis = drug_mice.index.values
y_axis = drug_mice.values

# Create a Pyplot bar plot based off of the group series from before and label the title
plt.bar(x_axis, y_axis, color='b', alpha=0.8, align='center')

# Set the xlabel and ylabel, title using class methods
plt.title("Number of Mice Tested per Treatment")
plt.xlabel("Drug Regimen")
plt.ylabel("Number of Mice")
plt.xticks(rotation=45)
plt.show()


# In[95]:


gender = clean_df["Sex"].value_counts()
gender


# In[ ]:





# In[96]:


# Generate a pie plot showing the distribution of female versus male mice using pyplot
labels = ['Female', 'Male']
plot = gender.plot.pie(y='Total Count', autopct="%1.0f%%")
plt.title('Mice Gender')
plt.ylabel('Sex')
plt.show()


# In[97]:


# Calculate the final tumor volume of each mouse across four of the treatment regimens:  
# Capomulin, Ramicane, Infubinol, and Ceftamin
Capomulin = clean_df.loc[clean_df["Drug Regimen"] == "Capomulin",:]
Infubinol = clean_df.loc[clean_df["Drug Regimen"] == "Infubinol", :]
Ceftamin = clean_df.loc[clean_df["Drug Regimen"] == "Ceftamin", :]
Ramicane = clean_df.loc[clean_df["Drug Regimen"] == "Ramicane", :]


# In[98]:


#Capomulin
Capomulin_last = Capomulin.groupby('Mouse ID').max()['Timepoint']
Capomulin_vol = pd.DataFrame(Capomulin_last)
Capomulin_merge = pd.merge(Capomulin_vol, clean_df, on=("Mouse ID","Timepoint"),how="left")
Capomulin_merge.head()


# In[99]:


#Calculate quartiles, IQR, lower bound, upper bound, and potential outliers for Capomulin
Capomulin_tumors = Capomulin_merge["Tumor Volume (mm3)"]

quartiles =Capomulin_tumors.quantile([.25,.5,.75])
lowerq = quartiles[0.25]
upperq = quartiles[0.75]
iqr = upperq-lowerq


print(f"The lower quartile of Capomulin tumors: {lowerq}")
print(f"The upper quartile of Capomulin tumors: {upperq}")
print(f"The interquartile range of Capomulin tumors: {iqr}")
print(f"The median of Capomulin tumors: {quartiles[0.5]} ")

lower_bound = lowerq - (1.5*iqr)
upper_bound = upperq + (1.5*iqr)

print(f"Values below {lower_bound} might be outliers.")
print(f"Values above {upper_bound} might be outliers.")


# In[100]:


Infubinol_last = Infubinol.groupby('Mouse ID').max()['Timepoint']
Infubinol_vol = pd.DataFrame(Infubinol_last)
Infubinol_merge = pd.merge(Infubinol_vol, clean_df, on=("Mouse ID","Timepoint"),how="left")
Infubinol_merge.head()


# In[101]:


#Calculate quartiles, IQR, lower bound, upper bound, and potential outliers for Infubino
Infubinol_tumors = Infubinol_merge["Tumor Volume (mm3)"]

quartiles =Infubinol_tumors.quantile([.25,.5,.75])
lowerq = quartiles[0.25]
upperq = quartiles[0.75]
iqr = upperq-lowerq


print(f"The lower quartile of Infubinol tumors is: {lowerq}")
print(f"The upper quartile of Infubinol tumors is: {upperq}")
print(f"The interquartile range of Infubinol tumors is: {iqr}")
print(f"The median of Infubinol tumors is: {quartiles[0.5]} ")

lower_bound = lowerq - (1.5*iqr)
upper_bound = upperq + (1.5*iqr)


print(f"Values below {lower_bound} might be outliers.")
print(f"Values above {upper_bound} might be outliers.")


# In[102]:


Ceftamin_last = Ceftamin.groupby('Mouse ID').max()['Timepoint']
Ceftamin_vol = pd.DataFrame(Ceftamin_last)
Ceftamin_merge = pd.merge(Ceftamin_vol, clean_df, on=("Mouse ID","Timepoint"),how="left")
Ceftamin_merge.head()


# In[103]:


#Calculate quartiles, IQR, lower bound, upper bound, and potential outliers for Ceftamin
Ceftamin_tumors = Ceftamin_merge["Tumor Volume (mm3)"]

quartiles = Ceftamin_tumors.quantile([.25,.5,.75])
lowerq = quartiles[0.25]
upperq = quartiles[0.75]
iqr = upperq-lowerq

print(f"The lower quartile of Ceftamin is: {lowerq}")
print(f"The upper quartile of Ceftamin is: {upperq}")
print(f"The interquartile range of Ceftamin is: {iqr}")
print(f"The the median of Ceftamin is: {quartiles[0.5]} ")

lower_bound = lowerq - (1.5*iqr)
upper_bound = upperq + (1.5*iqr)

print(f"Values below {lower_bound} could be outliers.")
print(f"Values above {upper_bound} could be outliers.")


# In[104]:


Ramicane_last = Ramicane.groupby('Mouse ID').max()['Timepoint']
Ramicane_vol = pd.DataFrame(Ramicane_last)
Ramicane_merge = pd.merge(Ramicane_vol, clean_df, on=("Mouse ID","Timepoint"),how="left")
Ramicane_merge.head()


# In[105]:


#Calculate quartiles, IQR, lower bound, upper bound, and potential outliers for Ramicane
Ramicane_tumors = Ramicane_merge["Tumor Volume (mm3)"]

quartiles =Ramicane_tumors.quantile([.25,.5,.75])
lowerq = quartiles[0.25]
upperq = quartiles[0.75]
iqr = upperq-lowerq


print(f"The lower quartile of Ramicane is: {lowerq}")
print(f"The upper quartile of Ramicane is: {upperq}")
print(f"The interquartile range of Ramicane is: {iqr}")
print(f"The median of Ramicane is: {quartiles[0.5]} ")

lower_bound = lowerq - (1.5*iqr)
upper_bound = upperq + (1.5*iqr)

print(f"Values below {lower_bound} might be outliers.")
print(f"Values above {upper_bound} might be outliers.")


# In[106]:


#Make a boxplot for all of the treatments
plot = [Capomulin_tumors, Ramicane_tumors, Infubinol_tumors, Ceftamin_tumors]
Regimen= ['Capomulin', 'Ramicane', 'Infubinol','Ceftamin']

fig1, ax1 = plt.subplots(figsize=(15, 10))
ax1.set_title('Mouse Tumor Volume',fontsize =20)
ax1.set_ylabel('End Tumor Volume (mm3)',fontsize = 15)
ax1.set_xlabel('Drug Regimen',fontsize = 15)
ax1.boxplot(plot, 0, 'gD', labels=Regimen, widths = 0.4, patch_artist=True,vert=True)

plt.ylim(0, 80)

plt.show()


# In[107]:


#Find a mouse that was treated with Capomulin
s185 = Capomulin.loc[Capomulin["Mouse ID"] == "s185",:]
s185.head()


# In[108]:


#Create line chart for a single mouse under Capomulin treatment
s185_line = s185["Timepoint"]
tumor = s185["Tumor Volume (mm3)"]

plt.title('Capomulin treatmeant of mouse s185')
plt.plot(s185_line, tumor,linewidth=3,)
plt.xlabel('Timepoint (Days)')
plt.ylabel('Tumor Volume (mm3)')

plt.show()


# In[109]:


#Create scatter plot for Capomulin
Cap_avg = Capomulin.groupby(['Mouse ID']).mean()
plt.scatter(Cap_avg['Weight (g)'],Cap_avg['Tumor Volume (mm3)'])
plt.xlabel('Weight (g)')
plt.ylabel('Average Tumor Volume (mm3)')

plt.show()


# In[110]:


# Calculate the correlation coefficient and linear regression model for mouse weight and average tumor volume for Capomulin
correlation=round(st.pearsonr(Cap_avg['Weight (g)'],Cap_avg['Tumor Volume (mm3)'])[0],2)
print(f"The correlation between mouse weight and average tumor volume is {correlation}")


# In[111]:


#Find the slope and intercept
stats=st.linregress(Cap_avg['Weight (g)'],Cap_avg['Tumor Volume (mm3)'])
stats


# In[112]:


#Set slope and intercept as variables
slope = 0.9544396890241045
intercept = 21.552160532685015


# In[113]:


#Plot the linear regression model on top of the previous line chart
y_values = Cap_avg['Weight (g)']*slope+intercept
plt.scatter(Cap_avg['Weight (g)'],Cap_avg['Tumor Volume (mm3)'])
plt.plot(Cap_avg['Weight (g)'],y_values,color="black")
plt.xlabel('Weight')
plt.ylabel('Average Tumor Volume')

plt.show()


# In[ ]:




