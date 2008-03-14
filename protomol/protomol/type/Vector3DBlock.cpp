#include <protomol/type/Vector3DBlock.h>
#include <protomol/type/Array.h>

using namespace std;
using namespace ProtoMol;
//____ Vector3DBlock

void Vector3DBlock::boundingbox(Vector3D &minbb, Vector3D &maxbb) const {
  if (size() <= 0) {
    minbb = Vector3D(Constant::MAXREAL, Constant::MAXREAL, Constant::MAXREAL);
    maxbb = Vector3D(-Constant::MAXREAL,
      -Constant::MAXREAL,
      -Constant::MAXREAL);
  }

  minbb = vec[0];
  maxbb = vec[0];

  const unsigned int count = size();
  for (unsigned int i = 1; i < count; ++i) {
    if (vec[i].x < minbb.x)
      minbb.x = vec[i].x;
    else if (vec[i].x > maxbb.x)
      maxbb.x = vec[i].x;
    if (vec[i].y < minbb.y)
      minbb.y = vec[i].y;
    else if (vec[i].y > maxbb.y)
      maxbb.y = vec[i].y;
    if (vec[i].z < minbb.z)
      minbb.z = vec[i].z;
    else if (vec[i].z > maxbb.z)
      maxbb.z = vec[i].z;
  }
}

Vector3D Vector3DBlock::sum() const {
  if (vec.empty())
    return Vector3D(0.0, 0.0, 0.0);

  //  Loop through all the atoms and remove the motion of center of mass
  //  using an advanced method -- Kahan's magic addition algorithm to get
  //  rid of round-off errors: Scientific Computing pp34.

  Vector3D sum(vec[0]);
  Vector3D tempC(0.0, 0.0, 0.0);
  for (unsigned int i = 1; i < size(); i++) {
    Vector3D tempX(vec[i]);
    Vector3D tempY(tempX - tempC);
    Vector3D tempT(sum + tempY);
    tempC = (tempT - sum) - tempY;
    sum = tempT;
  }

  return sum;
}

bool Vector3DBlock::fitplane(Vector3D &normal, Real &d, Real &err,
                             int limit) const {
  // Clear return values
  normal = Vector3D(0.0, 0.0, 0.0);
  d = 0.0;
  err = 0.0;

  // Oh no ...
  if (size() < 2)
    return false;

  // Special case, just a perfect plain
  if (size() == 3) {
    normal = (vec[0] - vec[1]).cross(vec[2] - vec[1]);
    d = -normal.dot(vec[0]);
  }
  // General case with more than three points ...
  else {
    int n = 4;
    int m = size();
    Array<Real, 2> v(ArraySizes((unsigned int)n + 1) ((unsigned int)n + 1));
    Array<Real, 2> a(ArraySizes((unsigned int)m + 1) ((unsigned int)n + 1));
    vector<Real> rv1(n + 1);
    vector<Real> w(n + 1);

    int flag, i, its, j, jj, k, l = 0, nm = 0;
    Real anorm, c, f, g, h, s, scale, x, y, z;
    for (unsigned int i = 0; i < size(); i++) {
      a[i + 1][1] = vec[i].x;
      a[i + 1][2] = vec[i].y;
      a[i + 1][3] = vec[i].z;
      a[i + 1][4] = 1.0;
    }

    // Singular Value Decomposition
    g = scale = anorm = 0.0;

    for (i = 1; i <= n; i++) {
      l = i + 1;
      rv1[i] = scale * g;
      g = s = scale = 0.0;
      if (i <= m) {
        for (k = i; k <= m; k++)
          scale += fabs(a[k][i]);

        if (scale) {
          for (k = i; k <= m; k++) {
            a[k][i] /= scale;
            s += a[k][i] * a[k][i];
          }

          f = a[i][i];
          g = -sign(sqrt(s), f);
          h = f * g - s;
          a[i][i] = f - g;
          for (j = l; j <= n; j++) {
            for (s = 0.0, k = i; k <= m; k++)
              s += a[k][i] * a[k][j];

            f = s / h;
            for (k = i; k <= m; k++)
              a[k][j] += f * a[k][i];
          }

          for (k = i; k <= m; k++)
            a[k][i] *= scale;
        }
      }
      w[i] = scale * g;
      g = s = scale = 0.0;
      if (i <= m && i != n) {
        for (k = l; k <= n; k++)
          scale += fabs(a[i][k]);

        if (scale) {
          for (k = l; k <= n; k++) {
            a[i][k] /= scale;
            s += a[i][k] * a[i][k];
          }

          f = a[i][l];
          g = -sign(sqrt(s), f);
          h = f * g - s;
          a[i][l] = f - g;
          for (k = l; k <= n; k++)
            rv1[k] = a[i][k] / h;

          for (j = l; j <= m; j++) {
            for (s = 0.0, k = l; k <= n; k++)
              s += a[j][k] * a[i][k];

            for (k = l; k <= n; k++)
              a[j][k] += s * rv1[k];
          }

          for (k = l; k <= n; k++)
            a[i][k] *= scale;
        }
      }
      anorm = max(anorm, (fabs(w[i]) + fabs(rv1[i])));
    }

    for (i = n; i >= 1; i--) {
      if (i < n) {
        if (g) {
          for (j = l; j <= n; j++)
            v[j][i] = (a[i][j] / a[i][l]) / g;

          for (j = l; j <= n; j++) {
            for (s = 0.0, k = l; k <= n; k++)
              s += a[i][k] * v[k][j];

            for (k = l; k <= n; k++)
              v[k][j] += s * v[k][i];
          }
        }
        for (j = l; j <= n; j++)
          v[i][j] = v[j][i] = 0.0;
      }
      v[i][i] = 1.0;
      g = rv1[i];
      l = i;
    }

    for (i = min(m, n); i >= 1; i--) {
      l = i + 1;
      g = w[i];
      for (j = l; j <= n; j++)
        a[i][j] = 0.0;

      if (g) {
        g = 1.0 / g;
        for (j = l; j <= n; j++) {
          for (s = 0.0, k = l; k <= m; k++)
            s += a[k][i] * a[k][j];

          f = (s / a[i][i]) * g;
          for (k = i; k <= m; k++)
            a[k][j] += f * a[k][i];
        }

        for (j = i; j <= m; j++)
          a[j][i] *= g;
      } else
        for (j = i; j <= m; j++)
          a[j][i] = 0.0;

      a[i][i] += 1.0;
    }

    for (k = n; k >= 1; k--)
      for (its = 1; its <= limit; its++) {
        flag = 1;
        for (l = k; l >= 1; l--) {
          nm = l - 1;
          if ((double)(fabs(rv1[l]) + anorm) == anorm) {
            flag = 0;
            break;
          }
          if ((double)(fabs(w[nm]) + anorm) == anorm)
            break;
        }

        if (flag) {
          c = 0.0;
          s = 1.0;
          for (i = l; i <= k; i++) {
            f = s * rv1[i];
            rv1[i] = c * rv1[i];
            if ((double)(fabs(f) + anorm) == anorm)
              break;
            g = w[i];
            h = norm(f, g);
            w[i] = h;
            h = 1.0 / h;
            c = g * h;
            s = -f * h;
            for (j = 1; j <= m; j++) {
              y = a[j][nm];
              z = a[j][i];
              a[j][nm] = y * c + z * s;
              a[j][i] = z * c - y * s;
            }
          }
        }
        z = w[k];
        if (l == k) {
          if (z < 0.0) {
            w[k] = -z;
            for (j = 1; j <= n; j++)
              v[j][k] = -v[j][k];
          }
          break;
        }
        if (its == limit)
          return false;
        x = w[l];
        nm = k - 1;
        y = w[nm];
        g = rv1[nm];
        h = rv1[k];
        f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
        g = norm(f, 1.0);
        f = ((x - z) * (x + z) + h * ((y / (f + sign(g, f))) - h)) / x;
        c = s = 1.0;
        for (j = l; j <= nm; j++) {
          i = j + 1;
          g = rv1[i];
          y = w[i];
          h = s * g;
          g = c * g;
          z = norm(f, h);
          rv1[j] = z;
          c = f / z;
          s = h / z;
          f = x * c + g * s;
          g = g * c - x * s;
          h = y * s;
          y *= c;
          for (jj = 1; jj <= n; jj++) {
            x = v[jj][j];
            z = v[jj][i];
            v[jj][j] = x * c + z * s;
            v[jj][i] = z * c - x * s;
          }

          z = norm(f, h);
          w[j] = z;
          if (z) {
            z = 1.0 / z;
            c = f * z;
            s = h * z;
          }
          f = c * g + s * y;
          x = c * y - s * g;
          for (jj = 1; jj <= m; jj++) {
            y = a[jj][j];
            z = a[jj][i];
            a[jj][j] = y * c + z * s;
            a[jj][i] = z * c - y * s;
          }
        }

        rv1[l] = 0.0;
        rv1[k] = f;
        w[k] = x;
      }

    // Find smallest singular value
    j = -1;
    for (i = 1; i <= n; i++)
      if (j < 0 || w[j] > w[i])
        j = i;

    // Get solution
    normal.x = v[1][j];
    normal.y = v[2][j];
    normal.z = v[3][j];
    d = v[4][j];
  }

  // Nomralize normal and plane constant, z >= 0.0
  Real nd = normal.norm();
  if (nd == 0.0) {
    normal = Vector3D(0.0, 0.0, 0.0);
    d = 0.0;
    err = 0.0;
    return false;
  }

  if (normal.z < 0.0) {
    normal = -normal;
    d = -d;
  }
  normal = normal / nd;
  d = d / nd;

  // The error
  for (unsigned int i = 0; i < size(); i++)
    err = power<2>(vec[i].dot(normal) + d);

  err = sqrt(err);

  return true;
}

