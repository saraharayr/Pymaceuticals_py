#!/usr/bin/env python
# coding: utf-8

# In[1]:


# Dependencies and Setup
import matplotlib.pyplot as plt
import pandas as pd
import scipy.stats as st
import numpy as np
from scipy.stats import linregress

# Study data files
mouse_metadata_path = "Desktop/Pymaceuticals_py/csv_data/Mouse_metadata.csv"
study_results_path = "Desktop/Pymaceuticals_py/csv_data/Study_results.csv"

# Read the mouse data and the study results
mouse_metadata = pd.read_csv(mouse_metadata_path)
study_results = pd.read_csv(study_results_path)

# Combine the data into a single dataset

full_dataset = pd.merge(mouse_metadata, study_results, how='outer', on='Mouse ID')


# Display the data table for preview
full_dataset.head()


# In[2]:


# Checking the number of mice.
mice = len(full_dataset["Mouse ID"].value_counts())
total_mice = pd.DataFrame({"Total Mice":[mice]})
total_mice.head()


# In[3]:


# Getting the duplicate mice by ID number that shows up for Mouse ID and Timepoint. 
dups_mice = full_dataset[full_dataset.duplicated(['Mouse ID', 'Timepoint'])]
round(dups_mice,2)


# In[4]:


# Create a clean DataFrame by dropping the duplicate mouse by its ID.

full_dataset.drop( full_dataset[full_dataset['Mouse ID']== 'g989'].index, inplace=True)
round(full_dataset,2)


# In[5]:


# Checking the number of mice in the clean DataFrame.
mice2 = len(full_dataset["Mouse ID"].value_counts())
clean_mice = pd.DataFrame({"No Duplicates Mice":[mice2]})
clean_mice.head()


# In[6]:


# Use groupby and summary statistical methods to calculate the following properties of each drug regimen: 
# mean, median, variance, standard deviation, and SEM of the tumor volume.

tumorv_mean = full_dataset.groupby("Drug Regimen").mean()["Tumor Volume (mm3)"]
tumorv_median = full_dataset.groupby("Drug Regimen").median()["Tumor Volume (mm3)"]
tumorv_var = full_dataset.groupby("Drug Regimen").var()["Tumor Volume (mm3)"]
tumorv_std = full_dataset.groupby("Drug Regimen").std()["Tumor Volume (mm3)"]
tumorv_sem = full_dataset.groupby("Drug Regimen").sem()["Tumor Volume (mm3)"]
                                                              
# Generate a summary statistics table of mean, median, variance, standard deviation, and SEM of the tumor volume for each regimen

summary_df = round(pd.DataFrame({"Tumor Volume Mean": tumorv_mean, "Tumor Volume Median": tumorv_median,
                                   "Tumor Volume Standard Deviation": tumorv_std, "Tumor Volume SEM": tumorv_sem}),2)

# Assemble the resulting series into a single summary dataframe.
summary_df


# In[7]:


# Generate a summary statistics table of mean, median, variance, standard deviation, and SEM of the tumor volume for each regimen
# Using the aggregation method, produce the same summary statistics in a single line

summary_line = round(full_dataset.groupby(['Drug Regimen'])['Tumor Volume (mm3)'].agg(["mean", "median", "var", "std", "sem"],axis="columns"),2)
summary_line 


# In[8]:


# Generate a bar plot showing the total number of timepoints for all mice tested for each drug regimen using Pandas.
timepoints = full_dataset.groupby(["Drug Regimen"])["Timepoint"].count().sort_values(ascending=False)
pd_barplot= timepoints.plot.bar(title="Mice per Drug Regimen Tested")
pd_barplot.set_xlabel("Drug Regimen")
pd_barplot.set_ylabel("Number of Mice Tested")
plt.show()


# In[9]:


# Generate a bar plot showing the total number of timepoints for all mice tested for each drug regimen using pyplot.
timepoints2 = full_dataset.groupby(["Drug Regimen"])["Timepoint"].count().sort_values(ascending=False)
pyplot_bar = pd.DataFrame(timepoints2)
pyplot_bar.plot.bar(legend=False)
plt.title("Mice per Drug Regimen Tested")
plt.xlabel("Drug Regimen")
plt.ylabel("Number of Mice Tested")
plt.show()


# In[10]:


# Generate a pie plot showing the distribution of female versus male mice using Pandas
pd_pie = full_dataset.groupby("Sex").nunique()["Mouse ID"]
pd_pie.plot.pie(title="Distribution of Female vs Male Mice")


# In[11]:


# Generate a pie plot showing the distribution of female versus male mice using pyplot
pyplot_pie = full_dataset.groupby("Sex").nunique()["Mouse ID"]
pyplot_piedf = pd.DataFrame(pyplot_pie)
plt.pie(pyplot_pie, labels=pyplot_pie.index)
plt.ylabel("Mouse ID")
plt.title("Distribution of Female vs Male Mice")
plt.show()


# In[12]:


# Calculate the final tumor volume of each mouse across four of the treatment regimens:  
# Capomulin, Ramicane, Infubinol, and Ceftamin

# Start by getting the last (greatest) timepoint for each mouse
last_timepoint = full_dataset.groupby("Mouse ID").max()["Timepoint"]
last_timepointdf = pd.DataFrame(last_timepoint)

# Merge this group df with the original dataframe to get the tumor volume at the last timepoint
merged_data = pd.merge(last_timepointdf,full_dataset, on=("Mouse ID","Timepoint"))
round(merged_data,2)


# In[13]:


# Put treatments into a list for for loop (and later for plot labels)
drugs = ["Capomulin", "Ramicane", "Infubinol", "Ceftamin"]

# Create empty list to fill with tumor vol data (for plotting)

vol_cap = []
vol_ram = []
vol_inf = []
vol_ceft = []

# Calculate the IQR and quantitatively determine if there are any potential outliers.   
# Locate the rows which contain mice on each drug and get the tumor volumes
   
for index, row in merged_data.iterrows():
    if row["Drug Regimen"] == "Capomulin":
        vol_cap.append(row["Tumor Volume (mm3)"])
    if row["Drug Regimen"] == drugs[1]:
        vol_ram.append(row["Tumor Volume (mm3)"])
    if row["Drug Regimen"] == drugs[2]:
        vol_inf.append(row["Tumor Volume (mm3)"])
    if row["Drug Regimen"] == drugs[3]:
        vol_ceft.append(row["Tumor Volume (mm3)"])
    
# add subset 
drugsdf = pd.DataFrame({"Capomulin": vol_cap,
                       drugs[1]: vol_ram,
                       drugs[2]: vol_inf,
                       drugs[3]: vol_ceft})   
    
# Determine outliers using upper and lower bounds

drugs_upper= drugsdf.max()
drug_lower = drugsdf.min()

round(drugsdf,3)


# In[14]:


#Calculate the quartiles and IRQ and quantitatively determine if there are any potential outliers across all four treatment regimens.

#Capomulin
quartiles_cap = drugsdf[drugs[0]].quantile([.25, .5, .75])
lowq_cap = quartiles_cap[.25]
medq_cap = quartiles_cap[.5]
upq_cap = quartiles_cap[.75]
IQR_cap = upq_cap-lowq_cap
lowb_cap = lowq_cap - (1.5*IQR_cap)
upb_cap = upq_cap + (1.5*IQR_cap)
max_cap = drugsdf[drugs[0]].max()
min_cap = drugsdf[drugs[0]].min()

#Ramicane
quartiles_ram = drugsdf["Ramicane"].quantile([.25, .5, .75])
lowq_ram = quartiles_ram[.25]
medq_ram = quartiles_ram[.5]
upq_ram = quartiles_ram[.75]
IQR_ram = upq_ram-lowq_ram
lowb_ram = lowq_ram - (1.5*IQR_ram)
upb_ram = upq_ram + (1.5*IQR_ram)
max_ram = drugsdf[drugs[1]].max()
min_ram = drugsdf[drugs[1]].min()

#Infubinol
quartiles_inf = drugsdf[drugs[2]].quantile([.25, .5, .75])
lowq_inf = quartiles_inf[.25]
medq_inf = quartiles_inf[.5]
upq_inf = quartiles_inf[.75]
IQR_inf = upq_inf-lowq_inf
lowb_inf = lowq_inf - (1.5*IQR_inf)
upb_inf = upq_inf + (1.5*IQR_inf)
max_inf = drugsdf[drugs[2]].max()
min_inf = drugsdf[drugs[2]].min()

#Ceftamin
quartiles_ceft = drugsdf[drugs[3]].quantile([.25, .5, .75])
lowq_ceft = quartiles_ceft[.25]
medq_ceft = quartiles_ceft[.5]
upq_ceft = quartiles_ceft[.75]
IQR_ceft = upq_ceft-lowq_ceft
lowb_ceft = lowq_ceft - (1.5*IQR_ceft)
upb_ceft = upq_ceft + (1.5*IQR_ceft)
max_ceft = drugsdf[drugs[3]].max()
min_ceft = drugsdf[drugs[3]].min()

outliers = pd.DataFrame({"Drugs": drugs,
                        "Lower Quartile":[lowq_cap, lowq_ram, lowq_inf, lowq_ceft],
                        "Upper Quartile":[upq_cap, upq_ram, upq_inf, upq_ceft],
                        "IQR":[IQR_cap, IQR_ram, IQR_inf, IQR_ceft],
                        "Median":[medq_cap, medq_ram, medq_inf, medq_ceft],
                        "Upper Bound": [upb_cap, upb_ram, upb_inf, upb_ceft],
                        "Lower Bound": [lowb_cap, lowb_ram, lowb_inf, lowb_ceft],
                        "Max": [max_cap, max_ram, max_inf, max_ceft],
                         "Min" : [min_cap, min_ram, min_inf, min_ceft]})
outliers


# In[15]:


# Generate a box plot of the final tumor volume of each mouse across four regimens of interest
boxplot_list = [drugsdf[drugs[0]],
                drugsdf[drugs[1]],
                drugsdf[drugs[2]],
                drugsdf[drugs[3]]]

fig1, ax = plt.subplots(figsize=(10,8))
ax.set_title("Tumor volume accross 4 regimes of interest") 
ax.set_xlabel("Treatments") 
ax.set_ylabel("Final Tumor Volume (mm3)") 
ax.boxplot(boxplot_list, 0, "gD")
plt.xticks([1,2,3,4], drugs) 
plt.show()


# In[16]:


# Generate a line plot of tumor volume vs. time point for a mouse treated with Capomulin

mouse = full_dataset.loc[full_dataset["Mouse ID"] == "b128"]

plt.figure(figsize=(10, 8), dpi=60)
plt.plot(mouse["Timepoint"], mouse["Tumor Volume (mm3)"], marker = "8")

plt.xlabel("Timepoint (days)")
plt.ylabel("Tumor Volume (mm3)")
plt.title("Capomulin Treatment of Mouse b128")
plt.show()


# In[17]:


# Generate a scatter plot of average tumor volume vs. mouse weight for the Capomulin regimen
capomulin_df = full_dataset.loc[full_dataset["Drug Regimen"] == "Capomulin"]

avg_m = full_dataset.groupby(["Mouse ID"]).mean()

plt.figure(figsize=(10, 8), dpi=60)
plt.scatter(avg_m["Weight (g)"],avg_m["Tumor Volume (mm3)"], s=50,c="blue")
plt.title("Average tumor volume vs. mouse weight for the Capomulin regimen")
plt.xlabel("Weight (g)")
plt.ylabel("Average Tumor Volume (mm3)")

plt.show()


# In[18]:


x_value = avg_m["Weight (g)"]
y_value = avg_m["Tumor Volume (mm3)"]
(slope, intercept, rvalue, pvalue, stderr) = linregress(x_value, y_value)
regress_values = x_value * slope + intercept

line_eq = f'y = {str(round(slope,2))}x + {str(round(intercept,2))}'
plt.scatter(x_value,y_value)
plt.plot(x_value,regress_values,"r-")
plt.annotate(line_eq,(26,35),fontsize=10,color="darkblue")
plt.title("Mouse weight vs. Avg. Tumor Volume")
plt.xlabel("Mouse weight (g)")
plt.ylabel("Tumor Volume (mm3)")
print(f"The r-squared is: {rvalue}")
print(f"The equation of the regression line is: {line_eq}")

plt.show()

