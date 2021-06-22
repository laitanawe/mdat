import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
def plot_fam(data,output):
    ax=sns.barplot(x=data['samples'],y=data['frequency'],hue=data['family'])
    ax.set(xlabel='Samples', ylabel='Frequency')
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., title='Families')
    plt.xticks(rotation=90)
    plt.show()
    ax.figure.savefig(f"{output}/families_barplot.png",dpi=1200,bbox_inches="tight")
    return ax
def plot_gen(data,output):
    ax=sns.barplot(x=data['samples'],y=data['frequency'],hue=data['genus'])
    ax.set(xlabel='Samples', ylabel='Frequency')
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., title='Genera')
    plt.xticks(rotation=90)
    plt.show()
    ax.figure.savefig(f"{output}/genera_barplot.png",dpi=1200,bbox_inches="tight")
    return ax
def plot_phyla(data,output):
    ax=sns.barplot(x=data['samples'],y=data['frequency'],hue=data['phylum'])
    ax.set(xlabel='Samples', ylabel='Frequency')
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., title='Phyla')
    plt.xticks(rotation=90)
    plt.show()
    ax.figure.savefig(f"{output}/phyla_barplot.png",dpi=1200,bbox_inches="tight")
    return ax
def plot_sp(data,output):
    ax=sns.barplot(x=data['samples'],y=data['frequency'],hue=data['species'])
    ax.set(xlabel='Samples', ylabel='Frequency')
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., title='Species')
    plt.xticks(rotation=90)
    plt.show()
    ax.figure.savefig(f"{output}/species_barplot.png",dpi=1200,bbox_inches="tight")
    return ax
def plot_strain(data,output):
    ax=sns.barplot(x=data['samples'],y=data['frequency'],hue=data['strain'])
    ax.set(xlabel='Samples', ylabel='Frequency')
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., title='Strains')
    plt.xticks(rotation=90)
    plt.show()
    ax.figure.savefig(f"{output}/strains_barplot.png",dpi=1200,bbox_inches="tight")
    return ax