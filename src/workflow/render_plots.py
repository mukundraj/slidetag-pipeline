import os


def get_tfmed_pts(tfm1, tfm2, ):


    # height=img_width_ss_tfmed # dimensions of transformed ss images exported from histolozee
    # width=img_height_ss_tfmed
    # # pts_tfmed_ss_nifti = [[-pt[0]*height, -pt[1]*width] for pt in pts_tfmed_ss_normalized]
    # pts_tfmed_ss_nifti = [[-pt[0]*width, -pt[1]*height] for pt in pts_tfmed_ss_normalized]


    # # warped_coords
    # nis_idx = str(nissl_id).zfill(3)
    # # cnissl_id = get_ofixed_nis_idx(nissl_id)

    # intrim_pts_file = tfm_folder+"/"+nis_idx+"_intrim_pts.txt"
    # with open(intrim_pts_file, 'w', newline='\n') as csvfile:
    #     writer = csv.writer(csvfile)
    #     for row in pts_tfmed_ss_nifti:
    #         writer.writerow(row)

    #     from_fiducials_file = tfm_folder+"/"+str(int(nis_idx))+"_f.csv"
    #     to_fiducials_file = tfm_folder+"/"+str(int(nis_idx))+"_t.csv"
    #     # nrrd_path = labelmap_folder+"/lmap_"+nis_idx+".nrrd"
    #     temp_output_csv_file = tfm_folder+"/"+nis_idx+"_tmp_output.csv"
    #     # call cpp script to transform output of previous step with fiduals/TPS, sample from nrrd and write output
    #     subprocess.run(["./build/cmapper2", 
    #                     from_fiducials_file,
    #                     to_fiducials_file,
    #                     intrim_pts_file,
    #                     "None",
    #                     nis_idx,
    #                     temp_output_csv_file])


    # ## read warped coords from temp file

    # warped_coords = []
    # with open(temp_output_csv_file, newline='\n') as csvfile:
    #     reader = csv.reader(csvfile)
    #     for row in reader:
    #         # print(row['first_name'], row['last_name'])
    #         row = list(map(float, row))
    #         point = [row[0], row[1]]
    #         warped_coords.append(point)

    # # frac_coords
    # width = img_width_nis_tfmed
    # height = img_height_nis_tfmed

    # frac_coords = [ [-float(pt[0])/width, -float(pt[1])/height] for pt in warped_coords ]
    pass


# transform coords based on T1 and T2

# read warped bead pos and transform to csv format

# plot tfmed points on nissl - keep bbox same as nissl dimensions

os.system(f'touch {snakemake.output.olayplot}')



