def my_legend(ax, loc='best', fancybox=False, framealpha=1, edgecolor='black',
              ncol=1, labelcolor='k', facecolor='w', fontsize=14, **kwargs):

    legend = ax.legend(loc=loc,
              fancybox=fancybox,
              framealpha=framealpha,
              edgecolor=edgecolor,
              fontsize=fontsize,
              ncol=ncol,
              labelcolor=labelcolor,
              facecolor=facecolor,
              **kwargs)

    legend.get_frame().set_linewidth(1.5)

    return legend



def my_text(ax, x, y, text, horizontalalignment='center',verticalalignment='center',
            bbox=dict(facecolor = 'white', edgecolor='black', pad=10.0), **kwargs):

    mytext = ax.text(x,
                    y,
                    text, 
                    horizontalalignment=horizontalalignment,
                    verticalalignment=verticalalignment,
                    transform=ax.transAxes,
                    bbox=bbox,
                    **kwargs)
        
    return mytext
