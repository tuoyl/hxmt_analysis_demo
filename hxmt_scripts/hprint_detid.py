#!/usr/bin/env python
 

print( "*********** Detector ID catalogue ***********" )

#payloads = raw_input("what do you want?(HE/ME/LE/ALL):\n")
payloads = 'ALL'
if payloads == "HE" or payloads == "he" or payloads == "ALL" or payloads == "all":
    print( "--------------------------")
    print( "--------------------------")
    print( "\nHE detector ID for SMALL FOV (5.7x1.1):")
    print( "0, 1, 3, 4, 5, 6, 7, 8, 10, 11, 12, 13, 14, 15, 17")
    print( "\nHE BLIND detector ID:")
    print( "16")
    print( "\nHE SMALL FOV information (5.7x1.1):")
    print( "FOV   0: 1, 3, 8, 11, 12")
    print( "FOV -60: 0, 4, 7,  9, 14, 17")
    print( "FOV +60: 2, 5, 6, 10, 13, 15")
    print( "\nHE LARGE FOV information (5.7x5.7):")
    print( "2, 9")
    print( "--------------------------")

   
if payloads == "ME" or payloads == "me" or payloads == "ALL" or payloads == "all":
    print( "--------------------------")
    print( "\nME ASIC ID for SMALL FOV (1x4):")
    print( "0-5, 7, 12-23, 25, 30-41, 43, 48-53")
    print( "\nME BLIND ASIC ID:")
    print( "10, 28, 46")
    print( "\nME ASIC ID with calibration source:")
    print( "6, 11, 24, 29, 42, 47")
    print( "\nME SMALL FOV information (1x4):")
    print( "FOV 0:   0-5,  7, 12-17")
    print( "FOV 1: 18-23, 25, 30-35")
    print( "FOV 2: 36-41, 43, 48-53")
    print( "\nME LARGE FOV information (4x4):")
    print( "8, 9, 26, 27, 44, 45")
    print( "--------------------------")

if payloads == "LE" or payloads == "le" or payloads == "ALL" or payloads == "all":
    print( "--------------------------")
    print( "\nLE detector ID for SMALL FOV (1.6x6):")
    print( '0,2-4,6-10,12,14,20,22-26,28,30,32,34-36,38-42,44,46,52,54-58,60-62,64,66-68,70-74,76,78,84,86,88-90,92-94')
    print( "\nLE BLIND detector ID for SMALL FOV:")
    print( "13, 45, 77")
    print( "\nLE BLIND detector ID for LARGE FOV (calibration source):")
    print( "21, 53, 85")
    print( "\nLE SMALL FOV information (1.6x6):")
    print( 'FOV 0:  0,  2-4, 6-10,12,14,20,22-26,28-30')
    print( 'FOV 1: 32,34-36,38-42,44,46,52,54-58,60-62')
    print( 'FOV 2: 64,66-68,70-74,76,78,84,86-90,92-94')
    print( "\nLE LARGE FOV information (6x6):")
    print( '1,5,11,15,27,31,33,37,43,47,59,63,65,69,75,79,91,95')
    print( "--------------------------")

