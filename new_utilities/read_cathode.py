#!/usr/bin/env python3

import shared_functions as sf

def main():

    usb_n = sf.get_usb_n(name="cathode")
    value = sf.read_hv(usb_n)
    print(f"Cathode at -{value} V")
    
    return 

if __name__ == "__main__":
    main()
