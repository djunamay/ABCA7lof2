import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
from statannotations.Annotator import Annotator

# functions
def compute_slope(point1, point2):
    m = (point2[1]-point1[1]) / (point2[0]-point1[0])
    return m
   
def get_points(last_point, data, measure):
    prev_point = last_point - 1
    point1 = np.array((np.array(data[data['Measurement']==prev_point]['Time']), np.array(data[data['Measurement']==prev_point][measure]))).reshape(-1)
    point2 = np.array((np.array(data[data['Measurement']==last_point]['Time']), np.array(data[data['Measurement']==last_point][measure]))).reshape(-1)
    return point1, point2

def compute_b(point, slope):
    b = point[1]-(slope*point[0])
    return b

def integrate(point1, point2, slope, b):
    v2 = (0.5*slope * point2[0]**2) + (b * point2[0])
    v1 = (0.5*slope * point1[0]**2) + (b * point1[0])
    return v2-v1

def compute_integral(temp, last_point, measure): 
    point1, point2 = get_points(last_point, temp, measure)
    slope = compute_slope(point1, point2)
    b = compute_b(point1, slope)
    return integrate(point1, point2, slope, b)

def compute_seahorse_measures_per_well_long_run(well_id, df_sub, value='OCR'):
    temp = df_sub[df_sub['Well']==well_id]
    integrals = [compute_integral(temp, x, value) for x in range(2,16)]
    
    
    non_mito = np.sum(integrals[12:14])
    basal = np.sum(integrals[0:2]) - non_mito
    etp_linked = np.sum(integrals[3:5]) - non_mito
    beta_abs = (basal-etp_linked)
    beta_ox = beta_abs/basal
    proton_leak = np.sum(integrals[6:8]) - non_mito
    atp_linked = basal-proton_leak
    max_resp = np.sum(integrals[9:11]) - non_mito
    
#     non_mito = np.sum(integrals[13:14])
#     basal = np.sum(integrals[1:2]) - non_mito
#     etp_linked = np.sum(integrals[4:5]) - non_mito
#     beta_abs = (basal-etp_linked)
#     beta_ox = beta_abs/basal
#     proton_leak = np.sum(integrals[7:8]) - non_mito
#     atp_linked = basal-proton_leak
#     max_resp = np.sum(integrals[10:11]) - non_mito
    
    
    CE = atp_linked/basal
    SRC = max_resp/basal
    ATP_of_MAX = atp_linked/max_resp

    return CE, SRC, ATP_of_MAX, basal, proton_leak, atp_linked, max_resp, beta_ox, etp_linked, beta_abs

def compute_seahorse_measures_per_well_short_run(well_id, df_sub):
    temp = df_sub[df_sub['Well']==well_id]
    integrals = [compute_integral(temp, x, 'OCR') for x in range(2,13)]
    non_mito = np.sum(integrals[9:11])
    basal = np.sum(integrals[0:2]) - non_mito
    proton_leak = np.sum(integrals[3:5]) - non_mito
    atp_linked = basal-proton_leak
    max_resp = np.sum(integrals[6:8]) - non_mito
    
    CE = atp_linked/basal
    SRC = max_resp/basal
    ATP_of_MAX = atp_linked/max_resp
#     temp = df_sub[df_sub['Well']==well_id]
#     non_mito = np.sum(temp[[x in set([10,11,12]) for x in temp['Measurement']]]['OCR']) #float(temp[temp['Measurement']==12]['OCR'])
#     basal = np.sum(temp[[x in set([1,2,3]) for x in temp['Measurement']]]['OCR'])- non_mito #float(temp[temp['Measurement']==3]['OCR'])
#     proton_leak = np.sum(temp[[x in set([4,5,6]) for x in temp['Measurement']]]['OCR'])- non_mito #float(temp[temp['Measurement']==6]['OCR']) 
#     atp_linked = basal-proton_leak
#     max_resp = np.sum(temp[[x in set([7,8,9]) for x in temp['Measurement']]]['OCR']) - non_mito #float(temp[temp['Measurement']==9]['OCR']) - non_mito

#     CE = atp_linked/basal
#     RCR = max_resp/proton_leak
#     ATP_of_MAX = atp_linked/max_resp

    return CE, SRC, ATP_of_MAX, basal, proton_leak, atp_linked, max_resp


def plot_boxplot_by_treatment(d, x_val, y_val, order, pairs, palette):
    
    ax = sns.boxplot(data = d, x = x_val, showfliers=False, y = y_val, palette = palette, dodge = True, order = order, width=.5, boxprops=dict(alpha=1), medianprops=dict(color='black', alpha=1), whiskerprops=dict(color='black', alpha=1), capprops=dict(color = 'black', alpha=1))
    sns.stripplot(data=d, x= x_val, y=y_val, dodge=True, jitter=True, alpha=1,  order = order, palette = palette)

     # Define pairs to compare
    annotator = Annotator(ax, pairs, data=d, x=x_val, y=y_val, order =order)
    annotator.configure(test='t-test_ind', text_format='full', loc='outside', verbose=2, show_test_name=False)
    
    annotator.apply_and_annotate()
    
    plt.legend(bbox_to_anchor=(1.04, 1), loc="upper left")
    plt.xticks(rotation=45)
    plt.xlabel('')
    
    
import os
from IPython.display import display,  clear_output
import ipywidgets as widgets
from PIL import Image

class ImageBrowser:
    def __init__(self, folder_path, long=False):
        self.folder_path = folder_path
        self.long = long
        
        self.df = pd.read_excel(self.folder_path, sheet_name='Rate')
        self.randomize_id()
          
        self.current_index = 0
        self.discard = []
        self.keep_ids = []
        
        self.previous = widgets.Button(description="previous")
        self.button = widgets.Button(description="next")
        
        self.warning = widgets.Button(description='discard', button_style='', value=False)
        self.warning.style.button_color = 'red' 
                
        self.keep = widgets.Button(description='keep', button_style='', value=False)
        self.keep.style.button_color = 'green' 
        
        self.output = widgets.Output()

        # Create horizontal boxes for the groups of buttons
        top_buttons = widgets.HBox([self.previous, self.button])
        bottom_buttons = widgets.HBox([self.warning, self.keep])

        # Stack the two horizontal boxes vertically
        button_layout = widgets.VBox([top_buttons, bottom_buttons])
        display(button_layout, self.output)


        self.previous.on_click(self.on_previous_clicked)
        self.button.on_click(self.on_button_clicked)
        self.warning.on_click(self.on_warning_clicked)
        self.keep.on_click(self.on_keep_clicked)

    def randomize_id(self):
        x = np.unique(self.df['Well'])
        randids = np.random.choice(np.arange(len(x)), size=len(x), replace=False)
        dictionary = dict(zip(x, randids))
        self.df['randID'] = [dictionary[x] for x in self.df['Well']]
        self.inv_dict = dict(zip(randids, x))
        
    
    def plot(self, ID, long):
        temp = self.df[self.df['randID']==ID]
        #print(temp)
        plt.figure()
        plt.plot(temp['Time'], temp['OCR'])
        plt.scatter(temp['Time'], temp['OCR'])

        if long:
            anno = ['Etomoxir', 'Olig', 'FCCP', 'Rot/Ant']
            for j,i in enumerate([2,5,8,11]):
                intervention = i
                y = np.max(temp['OCR'][temp['Measurement']==intervention])
                x = intervention*6.6 
                plt.annotate(
                anno[j],  # Text to display
                xy=(x, y),  # Point to annotate (tip of the arrow)
                xytext=(x, y + .1*y),  # Location of text (start of the arrow)
                arrowprops=dict(facecolor='red', arrowstyle='->'),  # Arrow properties
                horizontalalignment='center', verticalalignment='bottom' # Text alignment
                )
        else:
            anno = ['Olig', 'FCCP', 'Rot/Ant']
            for j,i in enumerate([2,5,8]):
                intervention = i
                y = np.max(temp['OCR'][temp['Measurement']==intervention])
                x = intervention*6.6 
                plt.annotate(
                anno[j],  # Text to display
                xy=(x, y),  # Point to annotate (tip of the arrow)
                xytext=(x, y + .1*y),  # Location of text (start of the arrow)
                arrowprops=dict(facecolor='red', arrowstyle='->'),  # Arrow properties
                horizontalalignment='center', verticalalignment='bottom' # Text alignment
                )

            #plt.close()
        return plt.gcf()
    
    def on_button_clicked(self,b):
        with self.output:
            clear_output(wait=True)
            
            display(self.plot(self.current_index, self.long))
            
            self.current_index +=1
            
    def on_previous_clicked(self,b):
        with self.output:
            clear_output(wait=True)

            self.current_index-=2
            
            self.discard = [i if i not in x.inv_dict[self.current_index] else None for i in self.discard]
            self.keep_ids = [i if i not in x.inv_dict[self.current_index] else None for i in self.keep_ids]

            display(self.plot(self.current_index))
            
    def on_warning_clicked(self,b):
        with self.output:
            (print('discarded'))
            self.discard.append(self.inv_dict[self.current_index-1])
            self.button.on_click(self.on_button_clicked)
    
    def on_keep_clicked(self,b):
        with self.output:
            (print('kept'))
            self.keep_ids.append(self.inv_dict[self.current_index-1])
