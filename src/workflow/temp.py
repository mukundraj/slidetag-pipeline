#! /usr/bin/env python

# https://simpleitk.org/SPIE2019_COURSE/01_spatial_transformations.html
# http://insightsoftwareconsortium.github.io/SimpleITK-Notebooks/Python_html/21_Transforms_and_Resampling.html
# https://medium.com/hipster-color-science/computing-2d-affine-transformations-using-only-matrix-multiplication-2ccb31b52181
# https://docs.scipy.org/doc/scipy/reference/generated/scipy.ndimage.affine_transform.html#scipy.ndimage.affine_transform
# https://stackoverflow.com/questions/39712767/how-to-set-size-for-scatter-plot

import SimpleITK as sitk
import numpy as np
from skimage.io import imread
from scipy import ndimage as ndi
from skimage.color import rgb2gray
from PIL import Image
import matplotlib.pyplot as plt

# reads ITK transform file and return extracted affine matrix
def read_tfm(tfm_path):
    tfm = sitk.ReadTransform(tfm_path)
    print(tfm)
    M = str(tfm).split('\n')[10:13]
    M = [i.strip().split(' ') for i in M]
    M = [[float(j) for j in i] for i in M]
    M = np.array(M)[:-1,:-1] # discarding 3rd dimension

    O = str(tfm).split('\n')[13][12:-1]
    O = [float(i) for i in O.split(',')]
    O = np.array(O[:-1]).reshape(2,1)

    C = str(tfm).split('\n')[14][12:-1]
    C = [float(i) for i in C.split(',')]
    C = np.array(C[:-1]).reshape(2,1)

    T = str(tfm).split('\n')[15][17:-1]
    T = [float(i) for i in T.split(',')]
    T = np.array(T[:-1]).reshape(2,1)

    A = np.zeros((3,3))
    A[:2,:2] = M
    A[:2,2] = T[:,0]
    return A


# transforms stag cell points, compute bbox and plot cropped image
def tfm_nissl(tfm_path, nissl_path, pts, op_path):
    A = read_tfm(tfm_path)
    # nissl_img = imread(nissl_path)
    # img = rgb2gray(imread(nissl_path))
    # A[0,2] = 1000
    # A[1,2] = 1000
    pts_tfmed = np.matmul(A, pts)

    # get bounding box of pts_tfmed
    bbox = np.zeros((2,2))
    bbox[0,0] = np.min(pts_tfmed[0,:])
    bbox[0,1] = np.max(pts_tfmed[0,:])
    bbox[1,0] = np.min(pts_tfmed[1,:])
    bbox[1,1] = np.max(pts_tfmed[1,:])
    bbox = bbox.astype(int)

    # print('pts_tfmed\n', pts_tfmed)
    # print('bbox\n', bbox)
    img = Image.open(nissl_path)
    area = (bbox[0,0], bbox[1,0], bbox[0,1], bbox[1,1])
    cropped_img = img.crop(area).convert('RGB')
    px = cropped_img.load()
    # px[0,0] = (0,0,0)
    # px[-1,1] = (255, 255, 255)
    print('channel\n', np.shape(cropped_img))
    # cropped_img.show()
    cropped_img.save(op_path, 'TIFF',dpi=(72,72))

# transforms cell positions and plots an image bounded by bbox of transformed cell positions
def tfm_stag(tfm_path, pts, op_path):
    A = read_tfm(tfm_path)
    pts_tfmed = np.matmul(A, pts)
    # get bounding box of pts_tfmed
    bbox = np.zeros((2,2))
    bbox[0,0] = np.min(pts_tfmed[0,:])
    bbox[0,1] = np.max(pts_tfmed[0,:])
    bbox[1,0] = np.min(pts_tfmed[1,:])
    bbox[1,1] = np.max(pts_tfmed[1,:])
    bbox = bbox.astype(int)
    # print('pts_tfmed2\n', pts_tfmed)
    # print('bbox2\n', bbox)

    # get bbox dimensions along both axes
    bbox_dims = np.zeros((2,1))
    bbox_dims[0,0] = bbox[0,1] - bbox[0,0]
    bbox_dims[1,0] = bbox[1,1] - bbox[1,0]

    print('bbox_dims\n', bbox_dims)

    # plot scatterplot of pts_tfmed with boundaries of bbox and no margins or padding of format tiff with no transparency and aspect ratio 1
    my_dpi = 72
    fig = plt.figure(figsize=(bbox_dims[0,0]/my_dpi, bbox_dims[1,0]/my_dpi), dpi=my_dpi, frameon=False)
    ax = plt.Axes(fig, [0., 0., 1, 1])
    ax.set_axis_off()
    fig.add_axes(ax)

    # plt.margins(0)
    # plt.axis('off')
    plt.scatter(pts_tfmed[0,:], pts_tfmed[1,:], s=500, c='#000000')
    plt.xlim(bbox[0,0], bbox[0,1])
    plt.ylim(bbox[1,0], bbox[1,1])

    plt.savefig(op_path, bbox_inches='tight', pad_inches=0, transparent=False, dpi='figure', format='tiff')
    im = Image.open(op_path)
    im.convert('RGB').save(op_path, 'TIFF',dpi=(72,72))



tfm_file ="/Users/mraj/Desktop/test/data/a5/rigid/nis_014_sl/Tfm1.txt"
# tfmed_pt = tfm.TransformPoint([0,0,0])
# print('A\n', A)

# create test 2D points in homogeneous coordinates
topr = 1050
test_pts = np.array([[0,0,1],[topr,0,1],[0,topr,1],[topr,topr,1]], dtype=np.float64).T
print('testpts\n', test_pts)

nis_file = "/Users/mraj/Desktop/test/data/a5/rigid/nis_014_sl/nis_014_sl.tif"
op_path = "/Users/mraj/Desktop/test/data/a5/rigid/nis_014_sl/nis_014_sl_crop.tif"
stag_file = ""
tfm_nissl(tfm_file, nis_file, test_pts, op_path)

op_path = "/Users/mraj/Desktop/test/data/a5/rigid/nis_014_sl/nis_014_sl_stag.tif"
tfm_stag(tfm_file, test_pts, op_path)

