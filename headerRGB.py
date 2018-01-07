import piexif
from PIL import Image
import piexif
import PIL.Image
from PIL import Image
from PIL.ExifTags import TAGS

def get_exif(fn):
    ret = {}
    i = Image.open(fn)
    info = i._getexif()
    
    for tag, value in info.items():
        decoded = TAGS.get(tag, tag)
        ret[decoded] = value
    return ret

filein = "science_29.jpg"
exif_dict = piexif.load(filein)
exif_dict['0th'][piexif.ImageIFD.DateTime] = '2017:04:23 22:30:30'
exif_bytes = piexif.dump(exif_dict)
piexif.insert(exif_bytes, filein)
