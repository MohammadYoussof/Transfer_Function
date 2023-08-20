import mod_tf as mtf

imgArr = mtf.Helper.LoadImageAsArray('simple.png')
res = mtf.MTF.CalculateMtf(imgArr, verbose=mtf.Verbosity.DETAIL)
