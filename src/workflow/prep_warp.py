"""
Script to geneate slicer proj files for warping

Created by Mukund on 2023-05-15
"""
from pathlib import Path
import sys
path_root = Path(__file__).parents[2]
sys.path.append(str(path_root))

import os
from src.python.slicer_mrml_gen import get_sub_text
import os
import numpy as np
import SimpleITK as sitk
import matplotlib.pyplot as plt
from PIL import Image
import PIL

# from produtils import dprint

# reads ITK transform file and return extracted affine matrix
# tfmed_yrange - if -1, indicates do operations on pts directly (needed initialy to determine tfmed_yrange)
#                if 0, indicates do operations after converting to slicer view space where [0,0] is at top left
def tfm_pts(tfm_path, pts, tfmed_yrange):

    tfm = sitk.ReadTransform(tfm_path)
    # print(tfm)
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

    # print(tfm.GetInverse().TransformPoint((0, 0, 0)))
    # print(tfm.GetInverse().TransformPoint((525, 0, 0)))
    # print(tfm.GetInverse().TransformPoint((525, 525, 0)))
    # print(tfm.GetInverse().TransformPoint((0, 525, 0)))
    # print(tfm.GetInverse().TransformPoint((341, 224, 0)))
    pts_tmp = pts[:2, :]

    # subtract y axis from 525 to get correct y axis
    if tfmed_yrange == -1:
        pts_tmp[1, :] =   pts_tmp[1, :]
    else:
        pts_tmp[1, :] =  tfmed_yrange - pts_tmp[1, :]

    pts_tfmed = []
    pts_tmp = np.vstack((pts_tmp, np.zeros((1, pts_tmp.shape[1]))))
    for i in range(pts_tmp.shape[1]):
        pts_tfmed.append(tfm.GetInverse().TransformPoint(pts_tmp[:, i].T))
    pts_tfmed = np.array(pts_tfmed).T


    # # print('M\n', M)
    # # print('O\n', O)
    # # print('C\n', C)
    # # print('T\n', T)
    # # print('A\n', A)

    # # print A
    # # print('A\n', A)
    # # invert A
    # A_inv = np.linalg.inv(A)
    # A_inv = np.array([[-2.09403934e-01, -8.25087080e-01,  5.27251472e+02],
    #                  [ 8.25087080e-01, -2.09403934e-01,  8.05165971e+02],
    #                  [ 0.00000000e+00,  0.00000000e+00,  1.00000000e+00]])
    # # A_inv = np.array([[-2.09403934e-01, -8.25087080e-01, 0],
    # #                  [ 8.25087080e-01, -2.09403934e-01, 0],
    # #                  [ 0.00000000e+00,  0.00000000e+00,  1.00000000e+00]])

    # # print('A_inv\n', A_inv)

    # pts_tfmed = np.matmul(A_inv, pts)
    return A, pts_tfmed

#  compute bbox of tfmed pts and plot cropped nissl image
def gen_and_save_cropped_nissl_img(nissl_path, op_path, bbox):

    # print('pts_tfmed\n', pts_tfmed)
    print('bbox\n', bbox)
    img = Image.open(nissl_path)
    area = (bbox[0,0], bbox[1,0], bbox[0,1], bbox[1,1])
    cropped_img = img.crop(area).convert('RGB')
    # cropped_img = cropped_img.transpose(PIL.Image.FLIP_LEFT_RIGHT)
    # cropped_img = cropped_img.transpose(PIL.Image.FLIP_TOP_BOTTOM)
    px = cropped_img.load()
    # px[0,0] = (0,0,0)
    # px[-1,1] = (255, 255, 255)
    print('channel\n', np.shape(cropped_img))
    # cropped_img.show()
    cropped_img.save(op_path, 'TIFF')
    # cropped_img.save(op_path, 'TIFF',dpi=(72,72))

# reads tfmed points and plots an image bounded by bbox of transformed cell positions
def gen_and_save_tfmed_stag_img(pts_tfmed, op_path, bbox):

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
    plt.scatter(pts_tfmed[0,:], pts_tfmed[1,:], s=5, c='#000000')
    # plt.scatter([526, 419, 849, 961], [807, 378, 263, 696], s=500, c=['red', 'green', 'blue', 'yellow'])
    plt.ylim(bbox[1,0], bbox[1,1])
    plt.xlim(bbox[0,0], bbox[0,1])
    ax = plt.gca()
    ax.set_ylim(ax.get_ylim()[::-1])
    # ax.set_xlim(ax.get_xlim()[::-1])
    print('ax.get_ylim()\n', ax.get_ylim())
    print('ax.get_xlim()\n', ax.get_xlim())

    plt.savefig(op_path, bbox_inches='tight', pad_inches=0, transparent=False, dpi='figure', format='tiff')
    im = Image.open(op_path)
    # im.convert('RGB').save(op_path, 'TIFF',dpi=(72,72))
    im.convert('RGB').save(op_path, 'TIFF')


def get_bbox(pts_tfmed):
    bbox = np.zeros((2,2))
    bbox[0,0] = np.min(pts_tfmed[0,:])# - 0.1*(np.max(pts_tfmed[0,:]) - np.min(pts_tfmed[0,:]))
    bbox[0,1] = np.max(pts_tfmed[0,:])# + 0.1*(np.max(pts_tfmed[0,:]) - np.min(pts_tfmed[0,:]))
    bbox[1,0] = np.min(pts_tfmed[1,:])# - 0.1*(np.max(pts_tfmed[1,:]) - np.min(pts_tfmed[1,:]))
    bbox[1,1] = np.max(pts_tfmed[1,:])# + 0.1*(np.max(pts_tfmed[1,:]) - np.min(pts_tfmed[1,:]))
    bbox = bbox.astype(int)
    return bbox


data_dir = snakemake.input.data
op_dir = snakemake.output.tfm1ed_pts_dir
os.system('mkdir -p ' + op_dir)
# read filenames in snakemake.input.data directory
filenames = os.listdir(snakemake.input.data)
# filter filenames to only include .csv files
filenames = [f for f in filenames if f.endswith('.csv')]

bbox = None
tfmed_yrange = None
# loop over filenames
for f in filenames:
    if f == 'bead_coords.csv':
        test_pts = np.genfromtxt(f'{data_dir}/{f}', delimiter=',', names=None, dtype=np.float64).T
        test_pts = test_pts[:2]

        # print test pts shape
        print('test_pts shape\n', test_pts.shape)
        # negate test points along y axis and y axis by subtracting from max value
        # test_pts[0,:] = np.max(test_pts[0,:]) - test_pts[0,:]
        # test_pts[1,:] = np.max(test_pts[1,:]) - test_pts[1,:]

        # swap x and y axes
        # test_pts = np.vstack((test_pts[1,:], test_pts[0,:]))
        
        # homogenize points
        test_pts = np.vstack((test_pts, np.ones((1, test_pts.shape[1]))))
        print('test_pts', test_pts)
        

        print('test_pts homo shape\n', test_pts.shape)
        
        bbbox_init = get_bbox(test_pts)
        print('bbbox_init', bbbox_init)

        # get bbox after tfm1
        A, tfmed_tst_pts = tfm_pts(snakemake.input.tfm1, test_pts, -1)
        tfmed_tst_pts_bbox = get_bbox(tfmed_tst_pts)
        tfmed_yrange = tfmed_tst_pts_bbox[1,1] - tfmed_tst_pts_bbox[1,0]
        print('tfmed_tst_pts bbox', tfmed_tst_pts_bbox)
        print('tfmed_tst_pts bbox y range', tfmed_yrange)
        # bbox_homo = np.array([[bbox]])


        # get tf1ed_bbox_dims

        # transform points and write to file
        A, pts_tfmed = tfm_pts(snakemake.input.tfm1, test_pts, tfmed_yrange)

        # print pts tfmed shape
        print('pts tfmed shape', np.shape(pts_tfmed[0,:]))

        print('min', np.min(pts_tfmed))
        print('max', np.argmax(pts_tfmed[0:]), np.max(pts_tfmed[0:]))
        print('max', np.argmax(pts_tfmed[1:]), np.max(pts_tfmed[1:]))

        bbox = get_bbox(pts_tfmed)
        print('bbbox\n', bbox)

        # save bbox as bbox.txt
        np.savetxt(snakemake.output.bbox, bbox, delimiter=",")

        # os.system('touch ' + op_file)
        # change file extension to .tif
        ip_nis_file = snakemake.input.nissl
        os.system('mkdir -p ' + snakemake.output.nis_imgs_dir)
        op_nis_name = snakemake.output.nis_imgs_dir+'/'+f.split('.')[0]+'_nissl.tif'
        print('op_nis_name\n', op_nis_name)
        gen_and_save_cropped_nissl_img(ip_nis_file, op_nis_name, bbox)

        # # subtract bbox[0,0] from all x coords and bbox[1,0] from all y coords
        # pts_tfmed[0,:] = pts_tfmed[0,:] - bbox[0,0]
        # pts_tfmed[1,:] = pts_tfmed[1,:] - bbox[1,0]

        # subtract y axis from 512  
        # pts_tfmed[1,:] = 512 - pts_tfmed[1,:]

        # pts_tfmed = -pts_tfmed

        # save tfmed points as .csv
        op_file_tfmed_pts = os.path.join(op_dir, f)
        np.savetxt(op_file_tfmed_pts, pts_tfmed[:2].T, delimiter=",")

# loop over filenames
# make stags imgs
os.system('mkdir -p ' + snakemake.output.stag_imgs_dir)
for f in filenames:
    # if f != 'bead_coords.csv':
    #     continue

    ip_file = os.path.join(data_dir, f)

    op_file = os.path.join(op_dir, f)

    # create tfmed points
    pts = np.genfromtxt(f'{ip_file}', delimiter=',', names=None, dtype=np.float64).T
    print(ip_file)
    print(np.shape(pts))
    test_pts = np.array([pts[0,:], pts[1,:], np.ones(np.shape(pts)[1])], dtype=np.float64)
    test_pts = test_pts[:2]
    # homogenize points
    test_pts = np.vstack((test_pts, np.ones((1, test_pts.shape[1]))))

    # transform points and write to file
    A, pts_tfmed = tfm_pts(snakemake.input.tfm1, test_pts, tfmed_yrange)

    # # subtract bbox[0,0] from all x coords and bbox[1,0] from all y coords
    # pts_tfmed[0,:] = pts_tfmed[0,:] - bbox[0,0]
    # pts_tfmed[1,:] = pts_tfmed[1,:] - bbox[1,0]

    # pts_tfmed = pts_tfmed[:2].T
    # # flip x and y axis by subtracting from max
    # pts_tfmed[0,:] = bbox[0,1] - pts_tfmed[0,:] + bbox[0,0]
    # pts_tfmed[1,:] = bbox[1,1] - pts_tfmed[1,:] + bbox[1,0]

    # # negate x and y axis
    # pts_tfmed[0,:] = -pts_tfmed[0,:]
    # pts_tfmed[1,:] = -pts_tfmed[1,:]
    # save tfmed points as .csv
    np.savetxt(op_file, pts_tfmed[:2].T, delimiter=",")


    op_stag_name = snakemake.output.stag_imgs_dir+'/'+f.split('.')[0]+'.tif'
    gen_and_save_tfmed_stag_img(pts_tfmed, op_stag_name, bbox)

# generate mrml file for warp from template

stag_imgs_dir = snakemake.output.stag_imgs_dir

# read filenames in stag_imgs_dir
stag_imgs = os.listdir(stag_imgs_dir)
stag_imgs = [f for f in stag_imgs if f.endswith('.tif')]

# remove extension from stag_imgs
stag_imgs = [f.split('.')[0] for f in stag_imgs]

# remove bead_coords from stag_imgs to avoid duplication    
stag_imgs = [f for f in stag_imgs if f != 'bead_coords']


sub_text = get_sub_text(stag_imgs)



# # copy template files
with open(snakemake.input.mrml_template, "rt") as fin:
    with open(f'{snakemake.output.mrml}', "wt") as fout:
        for line in fin:
            line = line.replace('fname_nissl', 'bead_coords_nissl')
            line = line.replace('fname_stag', 'bead_coords')
            line = line.replace('<Hidden></Hidden>', f'{sub_text}')
            fout.write(line)

os.system(f'cp {snakemake.input.from_fids} {snakemake.output.from_fids}')
os.system(f'cp {snakemake.input.to_fids} {snakemake.output.to_fids}')
# # os.system(f'touch {snakemake.output.tfm2}') # don't add this to output, empty file takes slicer to undesired state
