#!/usr/bin/env python3

start = 0
end   = 99
divisor=7
print("Printing out numbers from",start,"to",end, " not divisible by",divisor)


# list(range(0, 99+1))
for i in range(start, end+1): # need to add 1 to end to print "all the numbers from 0 to 99".  
    if i % divisor != 0:
        print(i)
    else: 
        continue
