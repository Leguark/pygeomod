'''
Created on Jul 1, 2015

@author: flow
'''

import numpy as np

import sys
import os
sys.path.append(os.path.realpath('/Users/Phisch/ownCloud/Philipp/Programming/Pygeomod/pygeomod-master/pygeomod'))

import geomodeller_xml_obj as GO
reload(GO)

class GeoSeis(GO.GeomodellerClass):
    '''
    classdocs
    '''


    def __init__(self):
        '''
        Constructor
        '''


    def get_points(self, **kwds):
        """Get Formation points from sections

        Points can either be defined with formation names or observation ids.
        In addition, one specific section can be defined (default: get points from all sections)


        **Optional Keywords**:
            - obs_id = string : observation id for points
            - formation_name = string : formation name
            - section_name = string : section name (default: use all sections)
            - return_obs_id = bool : return name of observation id
        """

        if kwds.has_key('debug') and kwds['debug']:
            debug = True
        else:
            debug = False

        try:
            self.section_dict
        except AttributeError:
            self.create_sections_dict()

        # list for original and new values, stored for all observation ids
        values = {}
        all_x_coords = []
        all_y_coords = []

        # outer loop: iterate through all sections
        for section in self.sections:
            sectionname = section.get("Name")
            # check if section name is specified as keyword:
            if kwds.has_key("section_name"):
                if sectionname != kwds['section_name']:
                    continue
        #    print("Extract points from section %s" % sectionname)
            # print("check section %s" % sectionname)
            points = self.get_formation_point_data(section)
            if points == None:
                print "No formation points defined in section " + sectionname
                continue # with the next section
            for point in points:
#                x_val = False
#                y_val = False
                obs_id = point.get("ObservationID")
#                data = point.find("{"+self.gml+"}Data")
                data = point.find("{"+self.xmlns+"}Data")
                formation_name = data.get("Name")
                fault_name = data.get("Name")
                # check if formation name is defined as keyword
                if kwds.has_key("formation_name") and (formation_name != kwds['formation_name']):
                    continue
                 # check if fault name is defined as keyword
                if kwds.has_key("fault_name") and (fault_name != kwds['fault_name']):
                    continue
                # Check observation id
                if kwds.has_key('obs_id') and (obs_id != kwds['obs_id']):
                    continue
#                print("Check obs_id %s" % obs_id)
                # add entry in values directory
                values[obs_id] = {}
                # get original point value
                element_point = point.find("{"+self.gml+"}LineString")
                element_coords = element_point.find("{"+self.gml+"}coordinates")
                point_list = element_coords.text.split(" ")
                if point_list[-1] == '':
                    point_list = point_list[0:-1]
                if len(point_list) > 1: # Multiple points defined by element
                    x_coords = []
                    y_coords = []
                    if debug:
                        print point_list
                    for point in point_list:
                        # if point == '': continue
                        a = point.split(',')
                        [x_coord, y_coord] = [float(a[0]), float(a[1])]
                        x_coords.append(x_coord)
                        y_coords.append(y_coord)
                        all_x_coords.append(x_coord)
                        all_y_coords.append(y_coord)
                    # convert to arrays for calculation
                    x_coord = np.array(x_coords)
                    y_coord = np.array(y_coords)
                obs_id_form = obs_id

        if kwds.has_key("return_obs_id") and kwds['return_obs_id']:
            return np.array(all_x_coords), np.array(all_y_coords), obs_id_form
        else:
            return np.array(all_x_coords), np.array(all_y_coords)

    def set_points(self, x_coords, y_coords, **kwds):
        """Set points in sections of Geomodeller model

        **Arguments**:
            - x_coords = list/ array : x-coordinate values
            - y_coords = list/ array : y-coordinate values

        **Optional keywords**:
            - obs_id = string : observation id for points
            - formation_name = string : formation name
            - section_name = string : section name (default: use all sections)

        """
        try:
            self.section_dict
        except AttributeError:
            self.create_sections_dict()

        # outer loop: iterate through all sections
        for section in self.sections:
            sectionname = section.get("Name")
            # check if section name is specified as keyword:
            if kwds.has_key("section_name"):
                if sectionname != kwds['section_name']:
                    continue
            #print("Extract points from section %s" % sectionname)
            # print("check section %s" % sectionname)
            points = self.get_formation_point_data(section)
            if points == None:
                print "No formation points defined in section " + sectionname
                continue # with the next section

            # set counter for array position
            array_count = 0

            for point in points:
                obs_id = point.get("ObservationID")
                element_point = point.find("{"+self.gml+"}LineString")

                data = point.find("{"+self.xmlns+"}Data")
                formation_name = data.get("Name")
                # check if formation name is defined as keyword
                if kwds.has_key("formation_name") and (formation_name != kwds['formation_name']):
                    continue

                # Check observation id
                if kwds.has_key('obs_id') and (obs_id != kwds['obs_id']):
                    continue

                element_coords = element_point.find("{"+self.gml+"}coordinates")
                point_list = element_coords.text.split(" ")
                if point_list[-1] == '':
                    point_list = point_list[0:-1]
                if len(point_list) > 1: # Multiple points defined by element
                    # now, reconstruct output format strings
                    out_text = ''
                    for _ in range(len(point_list)):


                        out_text += "%f,%f " % (x_coords[array_count], y_coords[array_count])
                        array_count += 1

                    element_coords.text = out_text

                else: # Single point only, use this to get started!
                    element_coords.text = "%f,%f" % (x_coords, y_coords)



if __name__ == '__main__':
    GS = GeoSeis()
    print GS
    GS.load_geomodeller_file("Pygeomodeller_Examples/models/Simple_Graben_3/Simple_graben_3.xml")
    x_coord, y_coord = GS.get_points(section_name = 'Section1')
    print(x_coord, y_coord)
    import matplotlib.pyplot as plt
    plt.plot(x_coord, y_coord, 'o')
    plt.show()
