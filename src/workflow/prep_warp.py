"""
Temporary scratch file for testing functions before adding to main pipeline

Created by Mukund on 2023-05-15
"""
import os
import numpy as np
import SimpleITK as sitk
import matplotlib.pyplot as plt
from PIL import Image

# from produtils import dprint

# reads ITK transform file and return extracted affine matrix
def tfm_pts(tfm_path, pts):
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
    A[2,2] = 1

    pts_tfmed = np.matmul(A, pts)
    return A, pts_tfmed

# transforms stag cell points, compute bbox and plot cropped image
def tfm_nissl(pts_tfmed, nissl_path, op_path):
    # nissl_img = imread(nissl_path)
    # img = rgb2gray(imread(nissl_path))
    # A[0,2] = 1000
    # A[1,2] = 1000
    # A, pts_tfmed = tfm_pts(tfm_path, pts)

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
def tfm_stag(pts_tfmed, op_path):
    # A, pts_tfmed = tfm_pts(tfm_path, pts)
    # pts_tfmed = np.matmul(A, pts)
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

# create test 2D points in homogeneous coordinates
topr = 1050
test_pts = np.array([[0,0,1],[topr,0,1],[0,topr,1],[topr,topr,1]], dtype=np.float64).T
print('testpts\n', test_pts)

print('tfm1', f'{snakemake.input.tfm1}')
# copy stag image after applying tfms and lift bbox
# os.system(f'touch {snakemake.output.nissl}')
# os.system(f'cp {snakemake.input.nissl} {snakemake.output.nissl}')

A, pts_tfmed = tfm_pts(snakemake.input.tfm1, test_pts)
print('pts_tfmed\n', pts_tfmed)
np.savetxt(snakemake.output.tfm1ed_pts, pts_tfmed[:2].T, delimiter=",")
tfm_nissl(pts_tfmed, snakemake.input.nissl, snakemake.output.nissl)

# copy nissl image after cropping bbox after transform stag coords
# os.system(f'touch {snakemake.output.stag}')
# os.system(f'cp {snakemake.input.stag} {snakemake.output.stag}')
tfm_stag(pts_tfmed, snakemake.output.stag)

# os.system(f'touch {snakemake.output.mrml}')
# os.system(f'touch {snakemake.output.from_fids}')
# os.system(f'touch {snakemake.output.to_fids}')

# copy template files
with open("./templates/warp/warp.mrml", "rt") as fin:
    with open(f'{snakemake.output.mrml}', "wt") as fout:
        for line in fin:
            line = line.replace('fname_nis', f'{snakemake.wildcards.fname}')
            line = line.replace('fname_stag', f'stag_{snakemake.wildcards.fname}')
            fout.write(line)

os.system(f'cp templates/warp/F.fcsv {snakemake.output.from_fids}')
os.system(f'cp templates/warp/T.fcsv {snakemake.output.to_fids}')
# os.system(f'touch {snakemake.output.tfm2}') # don't add this to output, empty file takes slicer to undesired state
