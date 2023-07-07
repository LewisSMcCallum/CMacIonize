#include <memory>

class BarnesHutTree {
public:
  BarnesHutTree(const CoordinateVector<double> &grid_dimensions,const CoordinateVector<double> box_anchor, double theta)
      : grid_dimensions(grid_dimensions),box_anchor(box_anchor),_theta(theta) {

        if (grid_dimensions[0] == grid_dimensions[1] && grid_dimensions[1] == grid_dimensions[2]) {
          roots.push_back(new TreeNode(grid_dimensions[0],grid_dimensions/2.0 + box_anchor));
        } else {
          double unit_block = grid_dimensions.min();
          uint_fast32_t need_blocksx = int(grid_dimensions[0]/unit_block);
          uint_fast32_t need_blocksy = int(grid_dimensions[1]/unit_block);
          uint_fast32_t need_blocksz = int(grid_dimensions[2]/unit_block);
          for (uint_fast32_t i=0;i<need_blocksx;i++) {
            for (uint_fast32_t j=0;j<need_blocksy;j++) {
              for (uint_fast32_t k=0;k<need_blocksz;k++) {
                CoordinateVector<double> new_mid((double)i*unit_block - unit_block/2.0,(double)j*unit_block - unit_block/2.0,(double)k*unit_block - unit_block/2.0);
                roots.push_back(new TreeNode(unit_block,new_mid + box_anchor));
              }
            }
          }


        }

  }


  void insert(const HydroDensitySubGrid::hydroiterator &cell) {
    CoordinateVector<double> cellpos = cell.get_cell_midpoint();
    for (auto& root : roots) {
      CoordinateVector<double> nodemax = root->node_midpoint;
      CoordinateVector<double> nodemin = root->node_midpoint;
      for (uint_fast32_t i=0;i<3;i++) {
        nodemax[i] += root->dimensions/2.0;
        nodemin[i] -= root->dimensions/2.0;
      }
      if (cellpos[0] < nodemax[0] && cellpos[0] > nodemin[0] && cellpos[1] < nodemax[1] && cellpos[1] > nodemin[1] &&
        cellpos[2] < nodemax[2] && cellpos[2] > nodemin[2]) {
          root->insert(cell);
          break;
        }
    }

  }

  CoordinateVector<double> compute_gravity(HydroDensitySubGrid::hydroiterator& cell) const {
    // Compute gravitational acceleration by traversing the octree

    CoordinateVector<double> acceleration(0.0,0.0,0.0);
    double starting_pot = cell.get_hydro_variables().get_gravitational_potential();
    for (auto& root : roots) {
      acceleration += root->compute_gravity_recursive(cell,_theta);

    }

    double end_pot = cell.get_hydro_variables().get_gravitational_potential();
    cell.get_hydro_variables().set_gravitational_potential(starting_pot + (end_pot-starting_pot)/2.0);


    return acceleration;



  }

  ~BarnesHutTree() {

    for (auto& root : roots) {
      root->delete_nodes();
    }



  }
  class TreeNode {
  public:
    double mass;
    CoordinateVector<double> center_of_mass;
    double dimensions;
    CoordinateVector<double> node_midpoint;
    std::unique_ptr<TreeNode> children[8];

    TreeNode(const double node_dimensions, const CoordinateVector<double> node_midpoint)
        : mass(0), center_of_mass(CoordinateVector<double>()), dimensions(node_dimensions), node_midpoint(node_midpoint) {
      for (int i = 0; i < 8; ++i) {
        children[i] = nullptr;
      }
    }

    void delete_nodes() {
      for (int i = 0; i < 8; ++i) {
        if (children[i]) {
          children[i]->delete_nodes();
          children[i].reset();

        }
      }
    }

    void insert(const HydroDensitySubGrid::hydroiterator &cell) {
      double resolution = std::pow(cell.get_volume(),1./3.);
      double cell_mass = cell.get_hydro_variables().get_conserved_mass();
      center_of_mass = (center_of_mass * mass + cell.get_cell_midpoint() * cell_mass) / (mass + cell_mass);
      mass += cell_mass;

      if (0.99*dimensions <= resolution) {
        //this is a single cell node (leaf node)
        // just add the mass
        return;

      } else if (children[0] == nullptr) {
        //cell hasnt been subdivided, divide it up, add the mass to this node then insert cell into correct child node
        subdivide();
        children[child_index(cell.get_cell_midpoint())]->insert(cell);
        return;

      } else {
        // child already exists, add mass here then insert cell into correct one.

        children[child_index(cell.get_cell_midpoint())]->insert(cell);
        return;

      }
    }

    int child_index(const CoordinateVector<double> &position) const {
      int index = 0;
      if (position.x() >= node_midpoint.x()) index |= 1;
      if (position.y() >= node_midpoint.y()) index |= 2;
      if (position.z() >= node_midpoint.z()) index |= 4;
      return index;
    }

    void subdivide() {
      for (int i = 0; i < 8; ++i) {
        double child_dimensions = dimensions / 2;
        CoordinateVector<double> child_midpoint = node_midpoint + CoordinateVector<double>((i & 1 ? 1 : -1) * child_dimensions / 2,
                                                                                           (i & 2 ? 1 : -1) * child_dimensions / 2,
                                                                                           (i & 4 ? 1 : -1) * child_dimensions / 2);
        children[i].reset(new TreeNode(child_dimensions, child_midpoint));
      }
    }


    CoordinateVector<double> compute_gravity_recursive(HydroDensitySubGrid::hydroiterator& cell, double theta) {
      CoordinateVector<double> position = cell.get_cell_midpoint();
      double distance = (position-node_midpoint).norm();


      if (dimensions/distance <= theta || children[0] == nullptr) {
        if (distance < 1.e-2) {
          //dont let cell interact with itself
          return 0.0;
        }
        double force = 6.6743e-11 * mass / std::pow(distance, 3);
        double pot = -6.6743e-11*mass*cell.get_hydro_variables().get_conserved_mass()/std::abs(distance);
        cell.get_hydro_variables().set_gravitational_potential(cell.get_hydro_variables().get_gravitational_potential() + pot);
        return (center_of_mass - position) * force;
      } else {
          CoordinateVector<double> acceleration(0, 0, 0);
          for (const auto& child : children) {
            if (child) {
              acceleration += child->compute_gravity_recursive(cell, theta);
            }
          }
          return acceleration;
      }
    }

  };

private:
  std::vector<TreeNode*> roots;
  CoordinateVector<double> grid_dimensions;
  CoordinateVector<double> box_anchor;
  double _theta;




};
