#!/usr/bin/env python3

import cv2
import numpy as np

def main_old(file1, file2, alpha_file, outfile):

    # Read 1st image
    img1 = cv2.imread(file1, cv2.IMREAD_UNCHANGED)
    background = img1.copy()

    # read image
    img = cv2.imread(file2, cv2.IMREAD_UNCHANGED)
    ht, wd = img.shape[:2]

    # extract the BGR component and the alpha component
    bgr_fg = img[:,:,0:3]
    if (img.shape[2] == 4):
        alpha_fg = img[:,:,3]
        print('fg contains alpha channel')
    else:
        alpha_fg = 255 * np.ones(img.shape[:2])
        print('fg does not contain alpha channel')

    if False:
        # compute pct% of ht and 100-pct% of ht
        pct = 50
        ht2 = int(ht*pct/100)
        ht3 = ht - ht2
        wd2 = int(wd * pct / 100)
        wd3 = wd - wd2

        # create opaque white image for left
        top = np.full((ht3,wd), 255, dtype=np.uint8)
        left = np.full((ht,wd3), 255, dtype=np.uint8)

        # create vertical gradient for bottom
        btm = np.linspace(255, 0, ht2, endpoint=True, dtype=np.uint8)
        btm = np.tile(btm, (wd,1))
        btm = np.transpose(btm)

        right = np.linspace(1.0, 0, wd2, endpoint=True)
        right = right * right * 255
        right = np.tile(right, (ht,1))

        # stack top and bottom
        alpha = np.hstack((left,right))

        # multiply alpha and alpha_fg
        alpha_new = (alpha * alpha_fg.astype(np.float64) / 255).clip(0,255).astype(np.uint8)

        # put alpha channel into image
        overlay = bgr_fg.copy()
        #overlay = cv2.cvtColor(overlay, cv2.COLOR_BGR2BGRA)
        #overlay[:,:,3] = alpha_new
        
        #background = cv2.cvtColor(background, cv2.COLOR_BGR2BGRA)
        #background[:,:,3] = 255 #1-alpha_new
        
        # Overlay?
        #result = cv2.addWeighted(background, 1, overlay, 1, 0)
        alpha_stack = np.stack([alpha_new, alpha_new, alpha_new], 2) / 255
        result = cv2.add((1-alpha_stack) * background, alpha_stack * bgr_fg)

    else:
        img_alpha = cv2.imread(alpha_file) #, cv2.IMREAD_UNCHANGED)
        alpha = img_alpha / 255.
        print(alpha.min(), alpha.max())
        result = cv2.add((1.0-alpha) * background, alpha * bgr_fg)


def main(file1, file2, alpha_file, outfile):
    background = cv2.imread(file1)[:,:,:3]
    foreground = cv2.imread(file2)[:,:,:3]
    alpha = cv2.imread(alpha_file)[:,:,:3]

    # Convert uint8 to float
    foreground = foreground.astype(float)
    background = background.astype(float)

    # Normalize the alpha mask to keep intensity between 0 and 1
    alpha = alpha.astype(float)/255

    # Multiply the foreground with the alpha matte
    foreground = cv2.multiply(alpha, foreground)
   
    # Multiply the background with ( 1 - alpha )
    background = cv2.multiply(1.0 - alpha, background)

    # Add the masked foreground and background.
    outImage = cv2.add(foreground, background)
    # save result
    cv2.imwrite(outfile, outImage)

if __name__ == '__main__':
    import plac; plac.call(main)
