#!/usr/bin/env python3

import shared_functions as sf

def main():

    usb_n = sf.get_usb_n(name="gate")
    value = sf.read_hv(usb_n)
    #print(f"Gate at -{value} V")
    print(value)

    return value 

if __name__ == "__main__":
    main()
