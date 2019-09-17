import myvariant

mv = myvariant.MyVariantInfo()

HGVS = []
filename = input('Please enter filename : ')
HGVS = list(myvariant.get_hgvs_from_vcf(filename))

index = 0

#Accept name input from user
name = input('Please enter an identifier for your file: ')
outputfile = open('HGVS_' + name + '.txt', 'w')

print('Writing to file...')
#Write to file
for ids in HGVS:
    outputfile.write(ids + '\n')
print('Done.')
    
outputfile.close()

