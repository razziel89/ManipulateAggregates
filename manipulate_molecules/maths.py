import numpy as np
import math as m

class MathsException(Exception):
    pass

class WrongDataFormatException(MathsException):
    pass

#rotation about an arbitrary axis going through the origin
def _rotation_matrix(axis,angle):
    result=np.zeros((3,3))
    u,v,w = axis/np.linalg.norm(axis)
    sin = m.sin(angle)
    cos = m.cos(angle)
    #source: http://inside.mines.edu/fs_home/gmurray/ArbitraryAxisRotation/
    #main diagonal
    result[0][0] = u*u + (v*v + w*w) * cos
    result[1][1] = v*v + (u*u + w*w) * cos
    result[2][2] = w*w + (u*u + v*v) * cos
    #upper triangle
    result[1][0] = u*v * (1-cos) + w*sin
    result[2][0] = u*w * (1-cos) - v*sin
    result[2][1] = v*w * (1-cos) + u*sin
    #lower triangle
    result[0][1] = u*v * (1-cos) - w*sin
    result[0][2] = u*w * (1-cos) + v*sin
    result[1][2] = v*w * (1-cos) - u*sin
    return result

def _dot(v1,v2):
    return np.sqrt(np.sum(v1*v2))

def _get_angle(v1,v2):
    n1=np.array(v1)
    n1/=np.linalg.norm(n1)
    n2=np.array(v2)
    n2/=np.linalg.norm(n2)
    return m.acos(_dot(n1,n2))

def _tensor_of_inertia(cloud):
    """
    Compute the tensor of inertia of a point cloud.

    cloud: the point cloud. Must be a numpy array of
           shape N,3
    """
    try:
        if cloud.shape[1]==3:
            result=np.zeros((3,3))
        else:
            raise WrongDataFormatException("Point cloud passed to _tensor_of_interia must consist of 3-element vectors.")
    except AttributeError as e:
        raise WrongDataFormatException("Point cloud passed to _tensor_of_interia must be numpy array (have the attribute shape)",e)
    print cloud[0]
    print np.linalg.norm(cloud,axis=1)[0]
    temp=np.sum(cloud*cloud,axis=0)
    result[0][0]=temp[1]+temp[2]
    result[1][1]=temp[0]+temp[2]
    result[2][2]=temp[0]+temp[1]
    temp=[-np.sum(cloud.T[0]*cloud.T[1]),-np.sum(cloud.T[0]*cloud.T[2]),-np.sum(cloud.T[1]*cloud.T[2])]
    result[0][1]=temp[0]
    result[1][0]=temp[0]
    result[0][2]=temp[1]
    result[2][0]=temp[1]
    result[1][2]=temp[2]
    result[2][1]=temp[2]
    return result

def align_clouds(reference,cloud,limit=None):
    """
    This function aligns all entries in the second argument to the entries
    in the first argument. Entries have to be three element vectors.
    Aligning means: both clouds have the same tensor of inertia and the
    same center (average of all vectors).

    Will return the adjusted cloud.

    reference: the point cloud that serves as a reference
    cloud: the cloud whose entries shall be adjusted
    limit: by default, all entries from cloud will be used
           to compute the current center and tensor of interia.
           If limit is not None, only this many elements 
           starting with the first will be used but the whole 
           list will be adjusted.
    """
    ref=np.array(reference)
    if not ref.shape[1]==3 or not len(ref.shape)==2:
        raise WrongDataFormatException("Shape of reference cloud must be N,3 (must consist of 3-element vectors).")
    c=np.array(cloud)
    if not c.shape[1]==3 or not len(ref.shape)==2:
        raise WrongDataFormatException("Shape of to-be-changed point-cloud must be M,3 (must consist of 3-element vectors).")
    if not isinstance(limit,int):
        raise ValueError("Given limit is no integer but it must be.")
    #get center of point clouds
    ref_center = np.mean(ref,axis=0)
    if limit==None:
        c_center   = np.mean(c  ,axis=0)
    else:
        c_center   = np.mean(c[:limit]  ,axis=0)
    #center both clouds around 0,0,0
    ref = ref-ref_center
    c   = c-c_center
    #get tensors of interia
    ref_inertia = _tensor_of_inertia(ref)
    if limit==None:
        c_inertia   = _tensor_of_inertia(c)
    else:
        c_inertia   = _tensor_of_inertia(c[:limit])
    tempvec = np.cross(c[0],ref[0])
    tempvec /= np.linalg.norm(tempvec)
    angle = _get_angle(c[0],ref[0])
    mat = _rotation_matrix(tempvec,angle)
    print
    print angle, tempvec, np.dot(mat,c[0]), ref[0], c[0]
    print
    return c-c_center+ref_center

