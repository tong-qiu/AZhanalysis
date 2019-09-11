import os
for filename in os.listdir("fetch"):
    os.rename("fetch/" + filename, "fetch/" + "hist-" + filename)
    print(filename)