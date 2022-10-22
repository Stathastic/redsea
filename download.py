import os
import gdown
date="20221021"

url="https://drive.google.com/drive/folders/1V5kMwzfgyg6F3jOW4DUAvPaSjQXahSkY?usp=sharing"
output="/home/jianj0c/dataset/redsea/{date}".format(date=date)

gdown.download_folder(url=url,output=output,use_cookies=False,quiet=False)

