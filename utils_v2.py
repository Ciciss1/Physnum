import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from matplotlib.lines import Line2D
import matplotlib.patches as mpatches
from matplotlib.collections import PathCollection


dark_blue = (0.0, 0.0, 0.4)
light_blue = (0.6, 0.8, 1.0)

dark_orange = (0.9, 0.5, 0.0)
light_orange = (1.0, 0.8, 0.2)

dark_green = (0, 0.3, 0.6)
light_green = (0.2, 0.7, 0.3)

dark_pink= (0.3,0.2,0.8)
light_pink = (1,0.5,0.9)

def convert_RGB_01(r,g,b) :
    r = r/255
    g = g/255
    b = b/255
    
    return (r,g,b)

blue = convert_RGB_01(59.0,117.0,175.0)
red = convert_RGB_01(197.0,58.0,50)
orange = convert_RGB_01(238.0,138.0,54.0)
purple = convert_RGB_01(118.0,45.0,121.0)
green = convert_RGB_01(145.0,188.0,134.0)
dark_yellow = convert_RGB_01(184.0,156.0,61.0)



def return_number_with_precision(number, n):
    if abs(number) < 10 ** (-n):
        formatted_number = 10**(-n)
    else:
        formatted_number = "{:.{}f}".format(number, n)
    return formatted_number

def linear_fit(x,y,color,
               
               #digits after coma
               precisions,
               
               #multiply a by 10^multiply
               multiply = 0,
               
               ax = plt, bounds = [0,0],
               
               a_err= 0, b_err = 0, linewidth = 3, plot = True) :
    x = remove_nan_values(x)
    y = remove_nan_values(y)
    n = len(x)
    
    abool = False
    bbool = False
    
    if a_err != 0 :
        abool = True
        
    if b_err != 0 :
        bbool = True
    
    a,b = np.polyfit(x,y,deg=1)

    if bounds != [0,0] : 
        x_fit = np.linspace(bounds[0],bounds[1],10)
    else : 
        x_fit = np.linspace(x.min(),x.max(),10)
        
    a_error = np.sqrt(np.sum((y-a*x-b)**2)/((n-2)*np.sum((x-np.mean(x))**2)))
    b_error = np.sqrt(np.sum(x**2)/n) * a_error
    
    y_fit = a*x_fit + b
    
    to_return = [a,a_error]
    
    
    if multiply != 0 : 
        a = a * pow(10,multiply)
        a_error = a_error * pow(10,multiply)
        b = b * pow(10,multiply)
        b_error = b_error * pow(10,multiply)
        
    if abool :
        a_error = a_err
        
    if bbool :
        b_error = b_err

    a_error = return_number_with_precision(a_error,precisions[0])
    b_error = return_number_with_precision(b_error,precisions[1])
    a = return_number_with_precision(a,precisions[0])
    b = return_number_with_precision(b,precisions[1])
    
    
    
    if multiply != 0 : 
        label = rf"$y = [({a} \pm {a_error})x + ({b} \pm {b_error})] \cdot 10^{-multiply}$"
    else : 
        label = rf"$y = ({a} \pm {a_error})x  + ({b} \pm {b_error}) $"
    
    if plot :
        ax.plot(x_fit,y_fit,color=color,label=label,linewidth = linewidth,linestyle="--")
    
    return to_return

def remove_nan_values(input_list):
    return np.array([value for value in input_list if not np.isnan(value)])

def scatter_multiple (colors, labels, X,Y,markers, hide = [], markersize = 10, ax = plt, alpha = 1) : 
    n = len(X)
    
    for i in range(n): 
        if not i in hide : 
            ax.scatter(X[i],Y[i],label = labels[i], marker = markers[i], color = colors[i], s=markersize, alpha=alpha) 

def errorbars_multiple (colors, labels, X,Y,markers, Y_error, X_error = [], ecolors = [], markersize = 20, capsize = 10, capthick = 1,ax = plt) : 
    n = len(X)
    
    if ecolors == [] : 
        ecolors = colors
        
    if X_error == [] : 
        for i in range(n): 
            ax.errorbar(X[i],Y[i],label = labels[i], yerr=Y_error[i], linestyle = '', marker=markers[i], capsize=capsize,capthick=capthick,ecolor=ecolors[i],color=colors[i],markersize=markersize)
    else : 
        for i in range(n):  
            ax.errorbar(X[i],Y[i],label = labels[i], xerr=X_error[i], yerr=Y_error[i], linestyle = '', marker=markers[i], capsize=capsize,capthick=capthick,ecolor=ecolors[i],color=colors[i],markersize=markersize)

def create_figure_and_apply_format(figsize,
                                   #label and axis settings
                                   xlabel, ylabel, xy_fontsize=22, tick_fontsize=18, 
                                   
                                   #grid or not
                                   grid_bool = True) : 
    
    fig = create_fig(figsize)
    
    ax = fig.gca()
    
    set_axis_and_tick_properties(ax,xlabel, ylabel, xy_fontsize, tick_fontsize)
    
    #set_legend_properties(ax, ncol, loc, fontsize, fontweight, fontstyle, text_color, border_color, border_linewidth)
    
    if grid_bool : 
        ax.grid()
    
    plt.tight_layout()
    
    return ax,fig
    

def create_fig(figsize) : 
    return plt.figure(figsize = figsize) 


def import_csv(filename = "CSV.csv") : 
    file_path = filename

    with open(file_path, 'r') as file:
        lines = file.readlines()

    modified_lines = [line.replace(',', '.') for line in lines]

    with open(file_path, 'w') as file:
        file.writelines(modified_lines)



    data = np.genfromtxt(file_path, delimiter=';', skip_header=1,  dtype=float)

    return data


def set_legend_properties(ax,colors,markers,markersize=15,ncol = 1, loc = "best", fontsize=13, fontweight='normal', fontstyle='italic', text_color='black', border_color='black', border_linewidth=2):
    handles, labels = ax.get_legend_handles_labels()
    
    i = 0
    j= 0
    
    proxy_artists = [0 for i in range(len(handles))]
    
    for h in handles :
        if type(h) == Line2D :
            proxy_artists[i] = h
        if type(h) ==  PathCollection:
            proxy_artists[i] = Line2D([0], [0], linestyle='', marker=markers[j], markersize=markersize, color=colors[j])
            j+=1
        i+=1
    
    
    legend = ax.legend(proxy_artists, labels, ncol = ncol, loc = loc)
    

    for label in legend.get_texts():
        label.set_fontsize(fontsize)
        label.set_fontweight(fontweight)
        label.set_fontstyle(fontstyle)
        label.set_color(text_color)

    legend.set_frame_on(True)
    legend.get_frame().set_linewidth(border_linewidth)
    legend.get_frame().set_edgecolor(border_color)
    
    plt.tight_layout()
 
def set_axis_and_tick_properties(ax,x_label, y_label, xy_fontsize, tick_fontsize):
    ax.set_xlabel(x_label, fontsize=xy_fontsize)
    ax.set_ylabel(y_label, fontsize=xy_fontsize)
    ax.tick_params('both', labelsize = tick_fontsize)


def save_png(fig, png_name = "test.png") : 
    fig.savefig("png/" + png_name)
    
def save_pdf(fig,pdf_name) : 
    fig.savefig("pdf/" + pdf_name)
    
def x_axis_divide(ax,divide = 1000) : 
    
    def divide_by(x, pos):
        return f'{x/divide:.1f}'
    formatter = FuncFormatter(divide_by)

    ax.xaxis.set_major_formatter(formatter)
        
def y_axis_divide(ax,divide = 1000) : 
    def divide_by(y,pos):
        return f'{y/divide:.1f}'
    formatter = FuncFormatter(divide_by)

    ax.yaxis.set_major_formatter(formatter)
    