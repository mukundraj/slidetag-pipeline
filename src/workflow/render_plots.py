import os
import numpy as np
import subprocess
import csv
import matplotlib.pyplot as plt


def get_tfmed_pts(tfm1ed_pts_file, from_fids_file, to_fids_file, bbox):
    print(f'{tfm1ed_pts_file}')
    print(f'{from_fids_file}')
    print(f'{to_fids_file}')


    # read tfm1ed_pts_file using np loadtxt
    tfm1ed_pts = np.loadtxt(tfm1ed_pts_file, delimiter=',')
    print('tfmed pts', tfm1ed_pts.shape)

    # get bounding box of tfm1ed_pts
    # bbox = np.array([np.min(tfm1ed_pts, axis=0), np.max(tfm1ed_pts, axis=0)])
    # print('bboxx', bbox)

    # get bbox dimensions
    bbox_dims = bbox[1] - bbox[0]
    print('bbox_dims', bbox_dims)
    width = bbox_dims[0]
    height = bbox_dims[1]


    # subtract bbox[0,0] from all x coords and bbox[1,0] from all y coords
    tfm1ed_pts[:,0] = tfm1ed_pts[:,0] - bbox[0,0]
    tfm1ed_pts[:,1] = tfm1ed_pts[:,1] - bbox[0,1]
    # tfm1ed_pts[:,0] = tfm1ed_pts[:,0] - 391.9#attn
    # tfm1ed_pts[:,1] = tfm1ed_pts[:,1] - 235.8
    print('minpts', np.min(tfm1ed_pts, axis=0))
    print('maxpts', np.max(tfm1ed_pts, axis=0))

    # invert y axis by subtracting from height and invert x axis by subtracting from width
    # tfm1ed_pts[:,1] = height - tfm1ed_pts[:,1]


    # height=img_width_ss_tfmed # dimensions of transformed ss images exported from histolozee
    # width=img_height_ss_tfmed
    # # pts_tfmed_ss_nifti = [[-pt[0]*height, -pt[1]*width] for pt in pts_tfmed_ss_normalized]
    # pts_tfmed_ss_nifti = [[-pt[0]*width, -pt[1]*height] for pt in tfm1ed_pts]]

    pts_tfmed_ss_nifti = np.zeros(tfm1ed_pts.shape)
    pts_tfmed_ss_nifti[:,0] = -1*tfm1ed_pts[:,0]
    pts_tfmed_ss_nifti[:,1] = -1*tfm1ed_pts[:,1]
    pts_tfmed_ss_nifti = list(pts_tfmed_ss_nifti)

    # warped_coords

    tfm_folder = os.path.split(tfm1ed_pts_file)[0]
    print('tfm_folder', tfm_folder)

    intrim_pts_file = tfm_folder+"/../intrim_pts.csv"
    with open(intrim_pts_file, 'w', newline='\n') as csvfile:
        writer = csv.writer(csvfile)
        for row in pts_tfmed_ss_nifti:
            writer.writerow(row)

    from_fiducials_file = from_fids_file
    to_fiducials_file = to_fids_file
    # nrrd_path = labelmap_folder+"/lmap_"+nis_idx+".nrrd"
    temp_output_csv_file = tfm_folder+"/../tmp_output.csv"
    # call cpp script to transform output of previous step with fiduals/TPS, sample from nrrd and write output
    # cmd = f'./build/cmapper2 {from_fiducials_file} {to_fiducials_file} {intrim_pts_file} None 0 {temp_output_csv_file}'
    # print(cmd)
    # os.system(cmd)
    subprocess.run(["./build/cmapper2", 
                    from_fiducials_file,
                    to_fiducials_file,
                    intrim_pts_file,
                    "None",
                    "0",
                    temp_output_csv_file])


    # ## read warped coords from temp file

    warped_coords = []
    with open(temp_output_csv_file, newline='\n') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            # print(row['first_name'], row['last_name'])
            row = list(map(float, row))
            point = [row[0], row[1]]
            warped_coords.append(point)

    # # frac_coords
    # width = img_width_nis_tfmed
    # height = img_height_nis_tfmed

    pos_warped_coords = [ [-float(pt[0]), -float(pt[1])] for pt in warped_coords ]
    pos_warped_coords = np.array(pos_warped_coords)

    return pos_warped_coords, bbox, bbox_dims, tfm1ed_pts

# read all files in tfm1ed_pts_dir
tfm1ed_pts_files = os.listdir(snakemake.input.tfm1ed_pts_dir)

# keep only csv files in tfm1ed_pts_files
tfm1ed_pts_files = [f for f in tfm1ed_pts_files if f.endswith('.csv')]

# remove .csv from tfm1ed_pts_files
tfm1ed_pts_files = [f.split('.')[0] for f in tfm1ed_pts_files]

# iterate over tfm1ed_pts_files
os.system(f'mkdir -p {snakemake.output.olay_plots_dir}')


for tfm1ed_pts_file in tfm1ed_pts_files:
    # if tfm1ed_pts_file != 'bead_coords':
    #     continue
    ip_fpath = os.path.join(snakemake.input.tfm1ed_pts_dir, tfm1ed_pts_file) + '.csv'
    op_fpath = os.path.join(snakemake.output.olay_plots_dir, tfm1ed_pts_file) + '.png'

    # remove '_coords' from op_fpath string
    op_fpath = op_fpath.replace('_coords', '')

    # os.system(f'touch {op_fpath}')

    # read bbox from bbox_file
    bbox_file = snakemake.input.bbox_file
    bbox = np.loadtxt(bbox_file, delimiter=',').T
    bbox_dims = bbox[1] - bbox[0]
    print('bbox_dims', bbox_dims)

    print(ip_fpath)
    pos_warped_coords, XXbbox, XXbbox_dims, tfm1ed_pts = get_tfmed_pts(ip_fpath, snakemake.input.from_fids, snakemake.input.to_fids, bbox)
    # pos_warped_coords = np.loadtxt(ip_fpath, delimiter=',')
    # print(pos_warped_coords)
    # pos_warped_coords = -1 * pos_warped_coords


    # invert y axis by subtracting from max y
    # pos_warped_coords[:,1] = bbox_dims[1] - pos_warped_coords[:,1]
    # pos_warped_coords[:,1] = bbox_dims[1] - pos_warped_coords[:,1]
    # invert x axis by subtracting from max x
    # pos_warped_coords[:,0] = bbox_dims[0] - pos_warped_coords[:,0]
    # print(pos_warped_coords)
    np.savetxt(op_fpath+'.txt', pos_warped_coords, delimiter=',')

    my_dpi = 72
    fig = plt.figure(figsize=(bbox_dims[0]/my_dpi, bbox_dims[1]/my_dpi), dpi=my_dpi, frameon=False)
    ax = plt.Axes(fig, [0., 0., 1, 1])
    ax.set_axis_off()
    fig.add_axes(ax)
    print('pos_warped_coords', pos_warped_coords)
    print('bbox', bbox)
    print('bbox_dims', bbox_dims)
    # bbox_dims = [988.1-391.9, 834.2-235.8] # attn
    # plt.scatter(pos_warped_coords[:,0], pos_warped_coords[:,1], s=500, c=v1, cmap='Greens')
    plt.scatter(pos_warped_coords[:,0], pos_warped_coords[:,1], s=50, c='red')
    # plt.scatter(tfm1ed_pts[:,0], tfm1ed_pts[:,1], s=50, c='red')
    # plt.xlim(bbox[0,0], bbox[1,0])
    # plt.ylim(bbox[0,1], bbox[1,1])
    plt.xlim(0, bbox_dims[0])
    plt.ylim(0, bbox_dims[1])
    ax = plt.gca()
    ax.set_ylim(ax.get_ylim()[::-1])

    tfmed_plot = f'{snakemake.output.olay_plots_dir}/tmp_tfmed_plot.png'
    plt.savefig(f'{tfmed_plot}', bbox_inches='tight', pad_inches=0, transparent=True, dpi='figure', format='png')
    plt.close()

    nissl_file = f'{snakemake.input.nissl_dir}/bead_coords_nissl.tif'
    cmd = f'convert {nissl_file} {tfmed_plot} -compose over -composite {op_fpath}'
    os.system(cmd)

exit(0)

# transform coords based on T2
pos_warped_coords, bbox, bbox_dims = get_tfmed_pts(snakemake.input.tfm1ed_pts, snakemake.input.from_fids, snakemake.input.to_fids)

# plot tfmed points - keep bbox same as nissl dimensions

# # plot pos_warped_coords with bounding box bbox
# pts = np.genfromtxt(f'{snakemake.input.data}', delimiter=',', names=True, dtype=np.float64).T
# v1 = (pts['v1'])
# print('v1', v1.shape)

op_folder = os.path.split(snakemake.input.tfm1ed_pts)[0]
my_dpi = 72
fig = plt.figure(figsize=(bbox_dims[0]/my_dpi, bbox_dims[1]/my_dpi), dpi=my_dpi, frameon=False)
ax = plt.Axes(fig, [0., 0., 1, 1])
ax.set_axis_off()
fig.add_axes(ax)
print('pos_warped_coords', pos_warped_coords)
# plt.scatter(pos_warped_coords[:,0], pos_warped_coords[:,1], s=500, c=v1, cmap='Greens')
plt.scatter(pos_warped_coords[:,0], pos_warped_coords[:,1], s=500,)
plt.xlim(bbox[0,0], bbox[0,1])
plt.ylim(bbox[1,0], bbox[1,1])

tfmed_plot = f'{op_folder}/tfmed_plt_plot.png'
plt.savefig(f'{tfmed_plot}', bbox_inches='tight', pad_inches=0, transparent=True, dpi='figure', format='png')


# create overlay image on nissl

# os.system(f'touch {snakemake.output.olayplot}')
# https://stackoverflow.com/questions/67577459/imagemagick-composite-two-images-and-fill-transparent-data-with-original-data
cmd = f'convert {snakemake.input.nissl} {tfmed_plot} -compose over -composite {snakemake.output.olayplot}'
os.system(cmd)


