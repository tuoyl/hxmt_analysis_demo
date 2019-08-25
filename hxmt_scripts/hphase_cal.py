#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division
from astropy.io import fits
import numpy as np 
import matplotlib.pyplot as plt
import sys
import argparse

try:
    # Python 2
    xrange
except NameError:
    # Python 3, xrange is now named range
    xrange = range

#read data file
def read_data(filename, colname='TDB'):
    hdulist = fits.open(filename)
    tdb = hdulist[1].data.field(colname)
    hdulist.close()
    return tdb

#read par file
def read_par(parname):
    pardata = open(parname,'r')
    stdpar = []
    parameters = np.zeros(13,dtype='longdouble')
    for par in pardata:
        par = par[0:(len(par)-1)]
        stdpar.append(par)
    pardata.close()
    for i in xrange(len(stdpar)):
        if stdpar[i][:6]=='PEPOCH':
            PEPOCH_lst = stdpar[i].split(' ');PEPOCH = [x for x in PEPOCH_lst if x is not ''][1]
            parameters[0] = np.longdouble(PEPOCH) 
        if stdpar[i][:2]=='F0': 
            F0_lst = stdpar[i].split(' ');F0 = [x for x in F0_lst if x is not ''][1]
            parameters[1] = np.longdouble(F0) 
        if stdpar[i][:2]=='F1':
            F1_lst = stdpar[i].split(' ');F1 = [x for x in F1_lst if x is not ''][1]
            parameters[2] = np.longdouble(F1)
        if stdpar[i][:2]=='F2':
            F2_lst = stdpar[i].split(' ');F2 = [x for x in F2_lst if x is not ''][1]
            parameters[3] = np.longdouble(F2)
        if stdpar[i][:2]=='F3':
            F3_lst = stdpar[i].split(' ');F3 = [x for x in F3_lst if x is not ''][1]
            parameters[4] = np.longdouble(F3)
        if stdpar[i][:2]=='F4':
            F4_lst = stdpar[i].split(' ');F4 = [x for x in F4_lst if x is not ''][1]
            parameters[5] = np.longdouble(F4)
        if stdpar[i][:2]=='F5':
            F5_lst = stdpar[i].split(' ');F5 = [x for x in F5_lst if x is not ''][1]
            parameters[6] = np.longdouble(F5)
        if stdpar[i][:2]=='F6':
            F6_lst = stdpar[i].split(' ');F6 = [x for x in F6_lst if x is not ''][1]
            parameters[7] = np.longdouble(F6)
        if stdpar[i][:2]=='F7':
            F7_lst = stdpar[i].split(' ');F7 = [x for x in F7_lst if x is not ''][1]
            parameters[8] = np.longdouble(F7)
        if stdpar[i][:2]=='F8':
            F8_lst = stdpar[i].split(' ');F8 = [x for x in F8_lst if x is not ''][1]
            parameters[9] = np.longdouble(F8)
        if stdpar[i][:2]=='F9':
            F9_lst = stdpar[i].split(' ');F9 = [x for x in F9_lst if x is not ''][1]
            parameters[10] = np.longdouble(F9)
        if stdpar[i][:5]=='START':
            START_lst = stdpar[i].split(' ');START = [x for x in START_lst if x is not ''][1]
            parameters[11] = np.longdouble(START) 
        if stdpar[i][:6]=='FINISH':
            FINISH_lst = stdpar[i].split(' ');FINISH = [x for x in FINISH_lst if x is not ''][1]
            parameters[12] = np.longdouble(FINISH) 

    print( "...finish reading ephemeris file...")
    return parameters

def check_subset(time, tstart, tfinish):
    if min(time) >= tstart and tfinish:
        pass
        #print("time is a subset of [START, FINISH]")
    else:
        print("WARNING: time set is not in parameter time range [START, FINISH]")

def phi_cal(time, parfile, instrument='hxmt'):
    if instrument == 'hxmt' or instrument == 'HXMT':
        MJDREFF = 0.0007660185
        MJDREFI = 55927
    elif instrument == 'fermi' or instrument == 'FERMI':
        MJDREFF = 0.00074287037037037
        MJDREFI = 51910

    #read parfile and parameters
    parameters = read_par(parfile)
    check_subset(time, parameters[11], parameters[12])
    PEPOCH = parameters[0]
    pepoch = (PEPOCH - MJDREFF - MJDREFI)*86400
    F0 = parameters[1]
    F1 = parameters[2]
    F2 = parameters[3]
    F3 = parameters[4]
    F4 = parameters[5]
    F5 = parameters[6]
    F6 = parameters[7]
    F7 = parameters[8]
    F8 = parameters[9]
    F9 = parameters[10]

    data = time
    t0 = pepoch # !!! set one reference point
    T0 = t0/86400 + MJDREFF + MJDREFI
    dt = t0 - pepoch 
    f0 = F0
    f1 = F1
    f2 = F2
    f3 = F3
    f4 = F4
    f5 = F5
    f6 = F6
    f7 = F7
    f8 = F8
    f9 = F9 
    #print 'periodic parameters',f0,f1,f2,f3,f4,f5,f6,f7,f8,f9

    phi = np.mod((data-t0)*f0 + (1/2)*((data-t0)**2)*f1 + (1/6)*((data-t0)**3)*f2 + (1/24)*((data-t0)**4)*f3 + (1/120)*((data-t0)**5)*f4 +
            (1/np.math.factorial(6))*((data-t0)**6)*f5 + (1/np.math.factorial(7))*((data-t0)**7)*f6 + (1/np.math.factorial(8))*((data-t0)**8)*f7 + 
            (1/np.math.factorial(9))*((data-t0)**9)*f8 + (1/np.math.factorial(10))*((data-t0)**10)*f9 ,1.0)
    #print phi
    print( "...processing..." )
    #print( phi )
    return phi

def cp_table(old_table):
    col_names = old_table.names
    col_type = old_table.formats
    cp_col = []
    for i in xrange(len(col_names)):
        cp_col.append( fits.Column(name=col_names[i], array=old_table.field(col_names[i]), format=col_type[i]) )
    new_table = fits.BinTableHDU.from_columns(cp_col)
    return new_table

def write_file(datafile, phi):
    #read old
    hdulist_old = fits.open(datafile)
    prim_hdr_old = hdulist_old[0].header
    prim_hdr_new = fits.PrimaryHDU(header=prim_hdr_old)
    hdr_all = []
    for i in xrange(len(hdulist_old)):
        hdr_all.append(hdulist_old[i].header)
    #modify main table
    table1 = hdulist_old[1].data
    col_names = table1.names
    col_type = table1.formats
    cp_col = []
    if 'Phase' not in col_names:
        print("...adding a column to event file...")
        new_col = fits.Column(name='Phase', array=phi, format='1D')
        for i in xrange(len(col_names)):
            cp_col.append( fits.Column(name=col_names[i], array=table1.field(col_names[i]), format=col_type[i]) )
        cp_col.append(new_col)
        new_table1 = fits.BinTableHDU.from_columns(cp_col)
        if len(hdulist_old) <=2:
            hdulist_new = fits.HDUList([prim_hdr_new, new_table1])
        else:
            table_rest = []
            for i in np.arange(2, len(hdulist_old), 1):
                table_rest.append(cp_table(hdulist_old[i].data))
            hdulist_new = fits.HDUList([prim_hdr_new, new_table1]+table_rest)
        #update header and info
        hdr_all[1]['TFIELDS'] = hdr_all[1]['TFIELDS'] + 1
        hdr_all[1]['NAXIS1'] = hdr_all[1]['NAXIS1'] + 8
        hdr_all[1].append('TTYPE' + str(hdr_all[1]['TFIELDS']), 'Phase', 'label for field')
        hdr_all[1].append('TFORM' + str(hdr_all[1]['TFIELDS']), '1D', 'format of field')
        for i in xrange(len(hdulist_new)):
            hdulist_new[i].header = hdr_all[i]
        hdulist_new.writeto(datafile, overwrite=True)
        hdulist_new.close()
        hdulist_old.close()
    else:
        print("...column Phase exists, overwrite the column...")
        table1['Phase'][:] = phi
        hdulist_old.writeto(datafile, overwrite=True)
        hdulist_old.close()
    print("...Success...")
    return 

def pass_argument():
    evtfile = []
    parfile = []
    colname = ['TDB']
    instrument = ['HXMT']
    len_arg = len(sys.argv)
    if len_arg <=1:
        raise IOError("Error input argument, RUN 'python cal_phase.py -h' for help")
    if sys.argv[1] == '-h' or sys.argv[1] == '--help':
        print( "#############" )
        print( '' )
        print( "EXAMPLE: python hphase_cal.py evtfile=eventfile.FITS parfile=ephemeris.par" )
        print( "" )
        print( "    evtfile: The Event file containing the column of the Barycenter corrected time" )
        print( "    parfile: The name of ephemeris file" )
        print( '    colname(default argument): The column name of Barycenter corrected time(the default value is "TDB"' )
        print( '    instrument(default argument): The name of Instrument(HXMT/FERMI) (the default value is "HXMT")' )
        print( '' )
        print( "############" )
        return False, False, False, False
    for i in xrange(len_arg):
        arg = sys.argv[i]
        if i == 0:continue
        if sys.argv[1] == '-h' or sys.argv[1] == '--help':
            pass
        else:
            arg_split = arg.split('=')
            argname = arg_split[0]
            argval  = arg_split[1]
            if argname not in ['evtfile', 'parfile', 'colname', 'instrument']:
                raise IOError("No such argument %s, RUN 'python cal_phase.py -h' for help"% argval)
            if argname == 'evtfile': 
                evtfile.append(argval)
            elif argname == 'parfile': 
                parfile.append(argval)
            elif argname == 'colname':
                colname[0] = argval
            elif argname == 'instrument':
                instrument[0] = argval
    if not evtfile or not parfile:
        raise IOError("Error input argument, RUN 'python cal_phase.py -h' for help")
        return False, False, False, False
    else:
        return evtfile[0], parfile[0], colname[0], instrument[0]

        

def main(evtfile, parfile, colname='TDB', instrument='HXMT'):
    time = read_data(evtfile,colname=colname)
    phi = phi_cal(time, parfile, instrument=instrument)
    write_file(evtfile, phi)
    return


if __name__ == "__main__":
    if len(sys.argv) <=1:
        print( "RUN 'python cal_phase.py --help' for help" )
    else:
        evtfile, parfile, colname, instrument = pass_argument()
        #print( evtfile, parfile, colname, instrument )
        if evtfile and parfile:
            main(evtfile, parfile, colname=colname, instrument=instrument)
