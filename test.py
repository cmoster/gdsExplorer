
import gdspy

class Polygon(gdspy.Polygon):
    def test(self):
        print 'new feture'


points =[[1,1],[1,2],[2,2],[2,1]]

p = Polygon(points)

print isinstance(p,Polygon)

p.test()

# print p