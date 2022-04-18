import numpy as np
import os
import vtk


def load_data_id(path):
    fopen = open(path, 'r', encoding='utf-8')
    lines = fopen.readlines()
    id_list = []
    linecount = 0

    for line in lines:
        if line == '\n':
            continue
        id_list.append(line.strip('\n'))
        linecount = linecount + 1
    fopen.close()
    return id_list


def check_and_create_dirs(dir_list):
    for dir in dir_list:
        if not os.path.exists(dir):
            os.makedirs(dir)
            print(dir + ' does not exist. Created.')


def save_off_points(points, path):
    with open(path, "w") as file:
        file.write("OFF\n")
        file.write(str(int(points.shape[0])) + " 0" + " 0\n")
        for i in range(points.shape[0]):
            file.write(
                str(float(points[i][0])) + " " + str(float(points[i][1])) + " " + str(float(points[i][2])) + "\n")


def save_off_mesh(v, f, path):
    with open(path, "w") as file:
        file.write("OFF\n")
        v_num = len(v)
        f_num = len(f)
        file.write(str(v_num) + " " + str(len(f)) + " " + str(0) + "\n")
        for j in range(v_num):
            file.write(str(float(v[j][0])) + " " + str(float(v[j][1])) + " " + str(float(v[j][2])) + "\n")
        for j in range(f_num):
            file.write("3 " + str(int(f[j][0])) + " " + str(int(f[j][1])) + " " + str(int(f[j][2])) + "\n")


def save_coff_points(points, colors, path):
    with open(path, "w") as file:
        file.write("COFF\n")
        file.write(str(int(points.shape[0])) + " 0" + " 0\n")
        for i in range(points.shape[0]):
            file.write(str(float(points[i][0])) + " " + str(float(points[i][1])) + " " + str(float(points[i][2])) + " ")
            file.write(str(colors[i][0]) + " " + str(colors[i][1]) + " " + str(colors[i][2]) + "\n")


def save_graph(v, A, path):
    with open(path, "w") as file:
        file.write("g line\n")
        v_num = len(v)
        for j in range(v_num):
            file.write("v " + str(float(v[j][0])) + " " + str(float(v[j][1])) + " " + str(float(v[j][2])) + "\n")
        file.write("g\n")

        # A is a symmetric matrix
        for j in range(v_num):
            for k in range(j + 1, v_num):
                if A[j][k] == 1:
                    file.write("l " + str(j + 1) + " " + str(k + 1) + "\n")

def rotationMatrix(A,B):
    def crossProduct(a,b):
        c = np.zeros((3,))
        c[0] = a[1]*b[2]-a[2]*b[1]
        c[1] = a[2]*b[0]-a[0]*b[2]
        c[2] = a[0]*b[1]-a[1]*b[0]
        c = c/(np.sqrt(c.dot(c.T)))
        return c

    rotationAxis = crossProduct(A,B)
    A = A/(np.sqrt(A.dot(A.T)))
    B = B/(np.sqrt(B.dot(B.T)))
    cosValue = A.dot(B.T)
    rotationAngle = np.arccos(cosValue)
    sinValue = np.sin(rotationAngle)

    matrix = np.zeros((3,3))
    matrix[0,0] = cosValue + rotationAxis[0]**2*(1-cosValue)
    matrix[0,1] = rotationAxis[0]*rotationAxis[1]*(1-cosValue) - rotationAxis[2]*sinValue
    matrix[0,2] = rotationAxis[0]*rotationAxis[2]*(1-cosValue) + rotationAxis[1]*sinValue
    matrix[1,0] = rotationAxis[0]*rotationAxis[1]*(1-cosValue) + rotationAxis[2]*sinValue
    matrix[1,1] = cosValue + rotationAxis[1]**2*(1-cosValue)
    matrix[1,2] = rotationAxis[1]*rotationAxis[2]*(1-cosValue) - rotationAxis[0]*sinValue
    matrix[2,0] = rotationAxis[0]*rotationAxis[2]*(1-cosValue) - rotationAxis[1]*sinValue
    matrix[2,1] = rotationAxis[1]*rotationAxis[2]*(1-cosValue) + rotationAxis[0]*sinValue
    matrix[2,2] = cosValue + rotationAxis[2]**2*(1-cosValue)
    return matrix

def save_final_obj(center,radius,faces,A,path):

    edges = []
    # print(len(center),len(A),center.shape)
    for j in range(len(center)):
        for k in range(j + 1, len(center)):
            if A[j][k] == 1:
                edges.append([j,k])

    ori_vector = np.array([0,0,1])
    sp_, sp_f = load_off('sphere16.off')
    vs = {}
    base = 1
    points = []
    faces_ = []
    with open(path, "w") as file:
        # file.write('v ' + str(0) + ' ' + str(0) + ' ' + str(0) + '\n')
        for i, edge in enumerate(edges):
            idx1,idx2 = edge[0],edge[1]
            v1,v2 = center[idx1],center[idx2]
            r1,r2 = radius[idx1],radius[idx2]
            dst_vector = v1 - v2
            dst_vector = dst_vector/(np.sqrt(dst_vector.dot(dst_vector.T)))
            matrix = rotationMatrix(ori_vector,dst_vector)
            sp_v = matrix.dot(sp_.T).T
            
            if idx1 not in vs:
                vs[idx1] = 0
                v1_ = sp_v[:65]*r1 + v1
                for m in range(v1_.shape[0]):
                    file.write('v ' + str(v1_[m][0]) + ' ' + str(v1_[m][1]) + ' ' + str(v1_[m][2]) + '\n')
                # base = m * sp_v.shape[0] + 1
                for j in range(49):
                    # file.write('f ' + str(sp_f[j][0] + base) + ' ' + str(sp_f[j][1] + base) + ' ' + str(sp_f[j][2] + base) + '\n')
                    if j==0:
                        for jj in range(15):
                            faces_.append([base,jj+base+1,jj+base+2])
                        faces_.append([base,16+base,base+1])
                    elif j%16==0:
                        faces_.append([j+base,j+16+base,j+base+1])
                        faces_.append([j+base,j+1+base,j+base-15])
                    else:
                        faces_.append([j+base,j+16+base,j+base+17])
                        faces_.append([j+base,j+17+base,j+base+1])
                base += 49
            else:
                v1_ = sp_v[49:65]*r1 + v1
                for m in range(v1_.shape[0]):
                    file.write('v ' + str(v1_[m][0]) + ' ' + str(v1_[m][1]) + ' ' + str(v1_[m][2]) + '\n')

            for j in range(15):
                faces_.append([j+base,j+16+base,j+base+17])
                faces_.append([j+base,j+17+base,j+base+1])
            j += 1
            faces_.append([j+base,j+16+base,j+base+1])
            faces_.append([j+base,j+1+base,j+base-15])
            base += 16

            if idx2 not in vs:
                v2_ = sp_v[49:]*r2 + v2
                vs[idx2] = 0
                for m in range(v2_.shape[0]):
                    file.write('v ' + str(v2_[m][0]) + ' ' + str(v2_[m][1]) + ' ' + str(v2_[m][2]) + '\n')
                for j in range(48):
                    # file.write('f ' + str(sp_f[j][0] + base) + ' ' + str(sp_f[j][1] + base) + ' ' + str(sp_f[j][2] + base) + '\n')
                    if j%16==15:
                        faces_.append([j+base,j+16+base,j+base+1])
                        faces_.append([j+base,j+1+base,j+base-15])
                        # pass
                    else:
                        faces_.append([j+base,j+16+base,j+base+17])
                        faces_.append([j+base,j+17+base,j+base+1])
                j=64
                for jj in range(15):
                    faces_.append([j+base,j+base-15+jj,j+base-16+jj])
                faces_.append([j+base,j+base-16,j+base-1])
                base += 49
            else:
                v2_ = sp_v[49:65]*r2 + v2
                for m in range(v2_.shape[0]):
                    file.write('v ' + str(v2_[m][0]) + ' ' + str(v2_[m][1]) + ' ' + str(v2_[m][2]) + '\n')
            base += 16

        for face in faces:
            p0,p1,p2 = center[face[0]],center[face[1]],center[face[2]]
            r0,r1,r2 = radius[face[0]],radius[face[1]],radius[face[2]]
            n = [0.,0.,0.]
            vtk.vtkTriangle().ComputeNormal(p0,p1,p2,n)
            n = np.array(n)
            points = [center[face[i]]+n*radius[face[i]] for i in range(3)]
            points += [center[face[i]]-n*radius[face[i]] for i in range(3)]
            for p in points:
                file.write('v ' + str(p[0]) + ' ' + str(p[1]) + ' ' + str(p[2]) + '\n')
            faces_.append([base+0,base+1,base+2])
            faces_.append([base+3,base+4,base+5])
            base += 6

        for face in faces_:
            file.write('f ' + str(face[0]) + ' ' + str(face[1]) + ' ' + str(face[2]) + '\n')



def save_spheres(center, radius, path):
    sp_v, sp_f = load_off('sphere16.off')

    with open(path, "w") as file:
        for i in range(center.shape[0]):
            v, r = center[i], radius[i]
            v_ = sp_v * r
            v_ = v_ + v
            for m in range(v_.shape[0]):
                file.write('v ' + str(v_[m][0]) + ' ' + str(v_[m][1]) + ' ' + str(v_[m][2]) + '\n')

        for m in range(center.shape[0]):
            base = m * sp_v.shape[0] + 1
            for j in range(sp_f.shape[0]):
                file.write(
                    'f ' + str(sp_f[j][0] + base) + ' ' + str(sp_f[j][1] + base) + ' ' + str(sp_f[j][2] + base) + '\n')


def save_skel_mesh(v, f, e, path_f, path_e):
    f_file = open(path_f, "w")
    e_file = open(path_e, "w")
    v_num = len(v)
    f_num = len(f)
    e_num = len(e)

    for j in range(v_num):
        f_file.write('v ' + str(float(v[j][0])) + " " + str(float(v[j][1])) + " " + str(float(v[j][2])) + "\n")
    for j in range(f_num):
        f_file.write("f " + str(int(f[j][0]) + 1) + " " + str(int(f[j][1]) + 1) + " " + str(int(f[j][2]) + 1) + "\n")

    for j in range(v_num):
        e_file.write('v ' + str(float(v[j][0])) + " " + str(float(v[j][1])) + " " + str(float(v[j][2])) + "\n")
    for j in range(e_num):
        e_file.write("l " + str(int(e[j][0]) + 1) + " " + str(int(e[j][1]) + 1) + "\n")

    f_file.close()
    e_file.close()


def save_skel_xyzr(v, r, path):
    file = open(path, "w")
    v_num = len(v)
    file.write(str(v_num) + "\n")
    for i in range(v_num):
        file.write(
            str(float(v[i][0])) + " " + str(float(v[i][1])) + " " + str(float(v[i][2])) + " " + str(float(r[i])) + "\n")
    file.close()


def save_colored_weights(path, shape_name, weights, samples):
    skel_num = weights.shape[0]
    sample_num = weights.shape[1]
    min_gray = 200
    for i in range(skel_num):
        colors = np.zeros((sample_num, 3)).astype(np.int)
        max_w = max(weights[i].tolist())
        for j in range(sample_num):
            color = min_gray - int((weights[i][j] / max_w) * min_gray)
            colors[j] = np.array([color, color, color], np.int)

        save_coff_points(samples, colors, path + str(shape_name) + '_' + str(i) + '_weight.off')


def load_off(path):
    fopen = open(path, 'r', encoding='utf-8')
    lines = fopen.readlines()
    linecount = 0
    pts = np.zeros((1, 3), np.float64)
    faces = np.zeros((1, 3), np.int)
    p_num = 0
    f_num = 0

    for line in lines:
        linecount = linecount + 1
        word = line.split()

        if linecount == 1:
            continue
        if linecount == 2:
            p_num = int(word[0])
            f_num = int(word[1])
            pts = np.zeros((p_num, 3), np.float)
            faces = np.zeros((f_num, 3), np.int)
        if linecount >= 3 and linecount < 3 + p_num:
            pts[linecount - 3, :] = np.float64(word[0:3])
        if linecount >= 3 + p_num:
            faces[linecount - 3 - p_num] = np.int32(word[1:4])

    fopen.close()
    return pts, faces


def load_ply_points(pc_filepath, expected_point=2000):
    fopen = open(pc_filepath, 'r', encoding='utf-8')
    lines = fopen.readlines()
    pts = np.zeros((expected_point, 3), np.float64)

    total_point = 0
    feed_point_count = 0

    start_point_data = False
    for line in lines:
        word = line.split()
        if word[0] == 'element' and word[1] == 'vertex':
            total_point = int(word[2])
            # if expected_point > total_point:
            #     pts = np.zeros((total_point, 3), np.float64)
            # continue

        if start_point_data == True:
            pts[feed_point_count, :] = np.float64(word[0:3])
            feed_point_count += 1

        if word[0] == 'end_header':
            start_point_data = True

        if feed_point_count >= expected_point:
            break

    fopen.close()
    return pts

if __name__ == '__main__':
    main()
