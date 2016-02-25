#!/usr/bin/env python
import os
import sys
import itertools
from subprocess import Popen, PIPE
import tempfile as tf

class ConfigError(Exception):
    pass

#defaul values for gnuplot
GNUPLOT_DEFAULT = {
                     "dash"   : "AUTO",
                     "color"  : "AUTO",
                     "lines"  : True,
                     "points" : False,
                     "title"  : "AUTO",
                     "xcol"   : 1,
                     "ycol"   : 2,
                     "zcol"   : 3,
                     "type"   : "UNKNOWN",
                     "head"   : False,
                     "linewidth" : 1.0,
                     "pointsize" : 1.0,
                     "solid"  : True,
                     "opacity": 1.0
                    }

class gnuplot():

    def __init__(self,filename,size=(20,12),linewidth=1,
            xrange=None,yrange=None,correct_mark=True,correct_mark_dist=0.001,
            fontsize=10,xlog=False,ylog=False,classic_colors=True,
            dimensions=2,font="Helvetica"):

        self.tempfiles = []
        self.GP = Popen(['gnuplot'], stdin=PIPE, stdout=sys.stderr, stderr=sys.stderr)
        self.f  = self.GP.stdin
        self.rectanglecount = 1
        self.dimensions = dimensions
        self.fontsize = fontsize
        self.font = font
        if self.dimensions not in (2,3):
            raise ConfigError("Wrong number of dimensions provided.")
        self.f.write("set term postscript eps enhanced colour \"%s,%d\" size %dcm,%dcm linewidth %d\n"%(font,fontsize,size[0],size[1],linewidth))
        if classic_colors:
            self.f.write("set colors classic\n")
        if xlog:
            self.f.write("set logscale x\n")
        if ylog:
            self.f.write("set logscale y\n")
        self.correct_mark = correct_mark
        if correct_mark:
            self.xmarks={}
            self.ymarks={}
            self.correct_mark_dist = correct_mark_dist
        self.xrange=xrange
        if self.xrange is not None:
            self.f.write("set xrange [%f:%f]\n"%(xrange[0],xrange[1]))
        self.yrange=yrange
        if self.yrange is not None:
            self.f.write("set yrange [%f:%f]\n"%(yrange[0],yrange[1]))
        self.f.write("set output \"%s.eps\"\n"%(filename))

    def _set_dict(self,dict):
        self.dict = dict
    
    def _get(self,key):
        if key in self.dict:
            return self.dict[key]
        elif key in GNUPLOT_DEFAULT:
            return GNUPLOT_DEFAULT[key]
        else:
            raise KeyError("Key %s not provided by the current dictionary and no default set."%(key))

    def set_title(self,title):
        self.f.write("set title \"%s\"\n"%(title))
    
    def emptyplot(self):
        if self.xrange is None or self.yrange is None:
            raise ConfigError("Cannot perform emtpy plot if either xrange or yrange is not set.")
        self.f.write("plot NaN notitle\n")
        self.postplot()

    def postplot(self):
        self.f.write("unset arrow\n")
        self.f.write("unset label\n")

    def unset(self,prop,oneprop=True):
        if oneprop:
            iterable = [prop]
        for p in iterable:
            self.f.write("unset %s\n"%(p))

    def set(self,prop,oneprop=True):
        if oneprop:
            iterable = [prop]
        for p in iterable:
            self.f.write("set %s\n"%(p))

    def lineplot(self,data):
        if isinstance(data,dict):
            dict_list = [data]
        else:
            try:
                if False in (isinstance(d,dict) for d in data):
                    raise TypeError("")
                else:
                    dict_list = data
            except TypeError:
                raise TypeError("Data for lineplot is neither a dictionary nor a list of dictionaries.")
        breakchar=", "
        if self.dimensions == 2:
            self.f.write("plot ")
        elif self.dimensions == 3:
            self.f.write("splot ")
        count = 1
        for f in dict_list:
            self._set_dict(f)
            if f == dict_list[-1]:
                breakchar = "\n"

            if not self._get("lines") and not self._get("points"):
                raise ValueError("At least one of 'lines' or 'points' has to be declared otherwise nothing would be plotted.")

            if self._get("type") == "function":
                self.f.write("%s "%(self._get("function")))
            elif self._get("type") == "file":
                self.f.write("\"%s\" u "%(self._get("file")))
                #x coloumn
                if isinstance(self._get("xcol"),int):
                    self.f.write("($%d):"%(self._get("xcol")))
                else:
                    self.f.write("(%s):"%(self._get("xcol")))
                #y coloumn
                if isinstance(self._get("ycol"),int):
                    self.f.write("($%d)"%(self._get("ycol")))
                else:
                    self.f.write("(%s)"%(self._get("ycol")))
                #z coloumn, if present
                if self.dimensions == 3:
                    if isinstance(self._get("zcol"),int):
                        self.f.write(":($%d) "%(self._get("zcol")))
                    else:
                        self.f.write(":(%s) "%(self._get("zcol")))
                self.f.write(" ")
            elif self._get("type") == "UNKNOWN":
                raise ValueError("No plot type provided. Missing key 'type' from dictionary.")
            else:
                raise ValueError("Unknown plot type: %s"%(f["type"]))

            if self._set_style(count,"line"):
                count += 1
            self.f.write(breakchar)
        if self.correct_mark:
            self.xmarks = {}
            self.ymarks = {}
        self.postplot()

    def data_to_file(self,data,formatstring=None,delete=True):
        tempfile = tf.NamedTemporaryFile(delete=False)
        if formatstring is None:
            for datum in data:
                tempfile.write("\t".join(map(str,datum))+"\n")
        else:
            for datum in data:
                tempfile.write(formatstring%tuple(datum))
        tempfile.close()
        self.tempfiles.append((delete,tempfile.name))
        return tempfile.name

    def _set_style(self,count,type):
        used_count = False
        style = ""
        if type == "line":
            style += "w "
            if self._get("lines"):
                style += "l"
            if self._get("points"):
                style += "p "
                style += "ps %f "%(self._get("pointsize"))
            style += " "
            dash = self._get("dash")
            if dash == "AUTO":
                style += "dt %d "%(count)
                used_count = True
            else:
                style += "dt %d "%(dash)
            color = self._get("color")
            if color == "AUTO":
                style += "lc %d "%(count)
                used_count = True
            else:
                style += "lc %d "%(color)
            width = self._get("linewidth")
            style += "lw %f "%(width)

            title = self._get("title")
            if title == "AUTO":
                pass
            elif title is None:
                style += "notitle "
            else:
                style += "title \"%s\" "%(title)
        elif type == "rectangle":
            color = self._get("color")
            if color == "AUTO" or color is None:
                style += "fc "
            else:
                style += "fc %s "%(color)
            style += "lt -1 lw 0 "
            if self._get("solid"):
                style += "fs solid %f "%(self._get("opacity"))
            if self._get("border") is None:
                style += "noborder "
        elif type == "arrow":
            dash = self._get("dash")
            if dash == "AUTO":
                style += "dt %d "%(count)
                used_count = True
            else:
                style += "dt %d "%(dash)
            color = self._get("color")
            if color == "AUTO":
                style += "lc %d "%(count)
                used_count = True
            else:
                style += "lc %d "%(color)
            width = self._get("linewidth")
            style += "lw %f "%(width)
            if not self._get("head"):
                style += "nohead "
        self.f.write(style)
        return used_count

    def set_xrange(self,start,stop):
        self.f.write("set xrange [%f:%f]\n"%(start,stop))
        self.xrange=(start,stop)

    def set_yrange(self,start,stop):
        self.f.write("set yrange [%f:%f]\n"%(start,stop))
        self.yrange=(start,stop)

    def set_zrange(self,start,stop):
        if dimensions == 3:
            self.f.write("set zrange [%f:%f]\n"%(start,stop))
            self.zrange=(start,stop)
        else:
            raise ConfigError("Cannot set zrange for non-3d plot.")

    def set_stick(self,pos,height,color,base=0.0,width=0.5):
        try:
            if len(pos) != len(height):
                print >>sys.stderr,"To print several sticks, the positions list and the height list must have the same number of elements."
            else:
                gen = ((p,h) for p,h in zip(pos,height))
        except TypeError:
            gen = [(pos,height)]
        for p,h in gen:
            self.mark(p,"x",width=width,color=color,rectangle=False,opacity=1.0,center=True,extent=(base,base+h))

    def set_background(self,opacities,colors=None,nr_fields=None,direction="x",extent=None):
        if direction == "x":
            samerange = self.xrange
            otherrange = self.yrange
        elif direction == "y":
            samerange = self.yrange
            otherrange = self.xrange
        else:
            raise ValueError("Unknown direction \"%s\", must be x or y."%(direction))
        if self.dimensions != 2:
            raise ConfigError("Cannot create background for non-2d plot.")
        if extent is None:
            if otherrange is None:
                raise ValueError("Cannot create background in %s-direction without other range being set."%(direction))
            else:
                extent = otherrange 
        if samerange is None:
            raise ValueError("Cannot create background in %s-direction without same range being set."%(direction))
        try:
            if colors is None:
                colors = [None] * len(opacities)
        except TypeError:
            opacities = [opacities]
            colors = [None]
        try:
            if len(colors) != len(opacities):
                raise ValueError("Cannot create background, colors and opacities do not have the same number of elements.")
            else:
                iterable = [(c,o) for c,o in zip(colors,opacities)]
        except TypeError:
            iterable = [(colours,opacities)]

        if nr_fields is None:
            nr_fields = len(colors)

        result = []
        width = 1.0*(samerange[1]-samerange[0])/(nr_fields)
        pos   = samerange[0]
        count = 0
        for color,opacity in itertools.cycle(iterable):
            self.mark(pos,direction,width=width,color=color,rectangle=True,center=False,opacity=opacity,extent=extent)
            result.append((pos,pos+width))
            pos += width
            count += 1
            if count == nr_fields:
                break
        return result

    def mark(self,pos,direction,width=0.5,color=None,rectangle=False,opacity=1.0,center=True,extent=None,label=None,zpos=None):
        if direction == "x":
            hererange = self.yrange
            heremarks = self.xmarks
            startpos  = lambda p,e: (p,e[0])
            endpos    = lambda p,e: (p,e[1])
        elif direction == "y":
            hererange = self.xrange
            heremarks = self.ymarks
            startpos  = lambda p,e: (e[0],p)
            endpos    = lambda p,e: (e[1],p)
        else:
            raise ValueError("Unknown direction \"%s\", must be x or y."%(direction))
        if self.dimensions != 2:
            if rectangle:
                raise ConfigError("Cannot set %smark as rectangle for non-2d plot."%(direction))
            elif self.dimensions == 3:
                if zpos is None:
                    raise ConfigError("Cannot set %smark as arrow for non-2d plot without zpos defined.")
            else:
                raise ConfigError("Fatal internal error: wrong number of dimensions set. However that happened.")
        if extent is None:
            if hererange is None:
                raise ValueError("Cannot create %smark without other range being set."%(direction))
            else:
                extent = hererange
        if not rectangle:
            if self.correct_mark:
                if pos in heremarks:
                    heremarks[pos] += 1
                    _pos = pos + self.correct_mark_dist*(heremarks[pos])
                else:
                    heremarks[pos] = 0
                    _pos = pos
            self.f.write("set arrow from %8.7E,%8.7E"%startpos(_pos,extent))
            if self.dimensions == 3:
                self.f.write(",%8.7E"%(zpos))
            self.f.write(" to %8.7E,%8.7E"%endpos(_pos,extent))
            if self.dimensions == 3:
                self.f.write(",%8.7E"%(zpos))
            self.f.write(" ")
            self._set_dict({"head":False,"color":color if color is not None else 0,"linewidth":width,"dash":1})
            self._set_style(heremarks[pos],"arrow")
            if opacity != 1.0:
                print >>sys.stderr,"WARNING: opacity unequal 100% set, but is ignored for xmark that is not a rectangle."
        else:
            if center:
                pos -= 0.5*width
            self.f.write("set obj %d rect from "%(self.rectanglecount))
            self.f.write("%f,%f to "%startpos(pos,extent))
            self.f.write("%f,%f "%endpos(pos+width,extent))
            self._set_dict({"color":color,"opacity":opacity,"border":None})
            self._set_style(self.rectanglecount,"rectangle")
            self.rectanglecount += 1
        self.f.write("\n")
        if label is not None:
            if isinstance(label,dict):
                label = [label]
            for l in label:
                if "where" in l:
                    where = l["where"]
                else:
                    where = "tl"
                if where == "tl":
                    labelpos = (startpos(pos,extent)[0]+l["offset"][0],startpos(pos,extent)[1]+l["offset"][1])
                    l["pivot"] = "left"
                elif where == "bl":
                    labelpos = (startpos(pos,extent)[0]+l["offset"][0],endpos(pos,extent)[1]+l["offset"][1])
                    l["pivot"] = "left"
                elif where == "tr":
                    labelpos = (endpos(pos,extent)[0]+l["offset"][0],startpos(pos,extent)[1]+l["offset"][1])
                    l["pivot"] = "right"
                elif where == "br":
                    labelpos = (endpos(pos,extent)[0]+l["offset"][0],endpos(pos,extent)[1]+l["offset"][1])
                    l["pivot"] = "right"
                elif where == "c":
                    labelpos = (0.5*(startpos(pos,extent)[0]+endpos(pos,extent)[0])+l["offset"][0],
                                0.5*(startpos(pos,extent)[1]+endpos(pos,extent)[1])+l["offset"][1])
                    l["pivot"] = "center"
                else:
                    raise ValueError("Wrong value for \"where\" in dictionary. Must be one of [\"tl\",\"bl\",\"tr\",\"br\",\"c\"] but is %s."%(where))
                l["position"] = labelpos
                self.set_label(l)

    def set_label(self,label):
        if "font" in label:
            font = label["font"]
        else:
            font = self.font
        if "fontsize" in label:
            fontsize = label["fontsize"]
        else:
            fontsize = self.fontsize
        if "pivot" in label:
            pivot = label["pivot"]
        else:
            pivot = "center"
        if "rotation" in label:
            rotation = label["rotation"]
        else:
            rotation = 0.0
        if "depth" in label:
            depth = label["depth"]
        else:
            depth = "front"
        self.f.write("set label \"%s\" at %f,%f font \"%s,%d\" %s rotate by %.2f %s\n"%(
                label["text"],label["position"][0],label["position"][1],font,fontsize,depth,rotation,pivot))

    def finalize(self,delete=True):
        if not self.f.closed:
            self.f.close()
            rc = self.GP.wait()
            if delete:
                for d,filename in self.tempfiles:
                    if d:
                        os.remove(filename)
            return rc
