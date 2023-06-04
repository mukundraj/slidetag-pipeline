"""
Functions to generate mrml file for 3D slicer using template

Created by Mukund on 2016-11-30
"""

print('bp1')
from bs4 import BeautifulSoup



vvd = """<VectorVolumeDisplay
  id="vtkMRMLVectorVolumeDisplayNode3" name="VectorVolumeDisplay_1" hideFromEditors="true" selectable="true" selected="false" color="0.9 0.9 0.3" edgeColor="0 0 0" selectedColor="1 0 0" selectedAmbient="0.4" ambient="0" diffuse="1" selectedSpecular="0.5" specular="0" power="1" opacity="1" sliceIntersectionOpacity="1" pointSize="1" lineWidth="1" representation="2" lighting="true" interpolation="1" shading="true" visibility="true" visibility2D="false" visibility3D="true" edgeVisibility="false" clipping="false" sliceIntersectionThickness="1" frontfaceCulling="false" backfaceCulling="false" scalarVisibility="false" vectorVisibility="false" tensorVisibility="false" interpolateTexture="false" scalarRangeFlag="UseData" scalarRange="0 100" colorNodeID="vtkMRMLColorTableNodeGrey" activeAttributeLocation="point" viewNodeRef="" folderDisplayOverrideAllowed="true" window="254" level="127" upperThreshold="32767" lowerThreshold="-32768" interpolate="1" windowLevelLocked="false" autoWindowLevel="1" applyThreshold="0" autoThreshold="0" visualizationMode="0" scalarMode="0" glyphMode="01" ></VectorVolumeDisplay>"""

vas = """ <VolumeArchetypeStorage
  id="vtkMRMLVolumeArchetypeStorageNode6" name="VolumeArchetypeStorage_1" hideFromEditors="true" selectable="true" selected="false" fileName="tifs/stags/stag_0.tif" fileListMember0="stag_0.tif" useCompression="1" defaultWriteFileExtension="nrrd" readState="0" writeState="4" centerImage="0" UseOrientationFromFile="1" ></VolumeArchetypeStorage>"""

vv = """ <VectorVolume
  id="vtkMRMLVectorVolumeNode3" name="stag_0" hideFromEditors="false" selectable="true" selected="false" references="display:vtkMRMLVectorVolumeDisplayNode3;storage:vtkMRMLVolumeArchetypeStorageNode6;" userTags="" spacing="1 1 1" origin="0 0 0" voxelVectorType="colorRGB" ijkToRASDirections="-1   0   0 0   -1   0 0 0 1 " measurementFrame="1   0   0 0   1   0 0 0 1 " order="-1" ></VectorVolume>"""


# xml = "<Hid></Hid>"
def get_sub_text(stag_imgs):

    all_img_str = ''
    for idx, img_name in enumerate(stag_imgs):

        soup1 = BeautifulSoup(vvd, 'xml')
        soup2 = BeautifulSoup(vas, 'xml')
        soup3 = BeautifulSoup(vv, 'xml')
        soup2.find('VolumeArchetypeStorage')['fileName'] = f'tifs/stags/{img_name}.tif'
        soup2.find('VolumeArchetypeStorage')['fileListMember0'] = f'{img_name}.tif'
        soup3.find('VectorVolume')['name'] = img_name

        # get id of soup1
        soup1_id = soup1.find('VectorVolumeDisplay')['id']
        soup2_id = soup2.find('VolumeArchetypeStorage')['id']

        # append idx to soup1_id and soup2_id
        soup1_id = soup1_id + '_'+ str(idx)
        soup2_id = soup2_id + '_'+ str(idx)

        # set id in soup1 and soup2
        soup1.find('VectorVolumeDisplay')['id'] = soup1_id
        soup2.find('VolumeArchetypeStorage')['id'] = soup2_id

        ref_str = f'references="display:{soup1_id};storage:{soup2_id};"'
        # set references in soup3
        soup3.find('VectorVolume')['references'] = ref_str

        soup1str = soup1.prettify().split('\n')[-2]
        soup2str = soup2.prettify().split('\n')[-2]
        soup3str = soup3.prettify().split('\n')[-2]
        # print('soup1str', soup1str)
        img_str = soup1str +'\n'+ soup2str +'\n' + soup3str


        # soup = BeautifulSoup(xml, 'xml')
        # soup.find('VolumeArchetypeStorage')['fileName'] = f'tifs/stags/{img_name}.tif'
        # soup.find('VolumeArchetypeStorage')['fileListMember0'] = f'{img_name}.tif'
        # print('soup', soup.find("VectorVolume"))
        # print('soup', soup)
        # soup.find('VectorVolume')['name'] = img_name

        # print(soup.prettify().split('\n')[-1])
        # print(str(soup).split("\n")[-1])
        # soup.find('VectorVolumeDisplay')['id'] = 'vtkMRMLVectorVolumeDisplayNode3'
        # print (soup.find('VectorVolumeDisplay')['id'])
        # print('stag_imgs', stag_imgs)
        # print('soup', soup)
        # img_str = xml.prettify()
        # img_str = xml
        # img_str = img_str.replace('fname_stag', img_name)
        # print(img_str)

        # append to all_img_str
        all_img_str += img_str
    # print('all_img_str', all_img_str)
    return all_img_str

def hello():
    print('hello')

