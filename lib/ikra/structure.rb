require "ikra/c_boolean"

module Ikra
  Fix        = Struct.new("Fix", :x, :y, :z)
  Volume     = Struct.new("Volume", :x, :y, :z)
  Periodic   = Struct.new("Periodic", :x, :y, :z)
  Coordinate = Struct.new("Coordinate", :x, :y, :z)
  MDLAtom    = Struct.new("MDLAtom", :type, :coordinate, :fix, :visible)
  XYZAtom    = Struct.new("XYZAtom", :type, :coordinate, :energy)

  def summation num
    num == 0 ? 0 : num + summation(num - 1)
  end

  class MDL
    attr_accessor :condition, :volume, :size, :atoms

    def read fp
      @condition     = Condition.new
      @volume, @size = @condition.read_from_mdl(fp)
      @atoms         = []
      fp.each do |line|
        words           = line.split
        atom            = MDLAtom.new
        atom.type       = words[0]
        atom.coordinate = Coordinate.new(words[1].to_f, words[2].to_f, words[3].to_f)
        atom.fix        = Fix.new(words[4].to_i, words[5].to_i, words[6].to_i)
        atom.visible    = words[7].to_i
        @atoms.push atom
      end
    end

    def write fp
      @condition.write_to_mdl fp, @volume, @size
      @atoms.each do |atom|
        fp.puts "#{atom.type} % 20.15e % 20.15e % 20.15e #{atom.fix.x} #{atom.fix.y} #{atom.fix.z} #{atom.visible}"%[atom.coordinate.x, atom.coordinate.y, atom.coordinate.z]
      end
    end

    def replace_condition c
      # non replaced condition
      c.displacement_ux = Marshal.load(Marshal.dump(@condition.displacement_ux))
      c.displacement_uz = Marshal.load(Marshal.dump(@condition.displacement_uz))
      c.height          = Marshal.load(Marshal.dump(@condition.height))
      c.stress_ex       = Marshal.load(Marshal.dump(@condition.stress_ex))
      c.stress_ez       = Marshal.load(Marshal.dump(@condition.stress_ez))
      c.spbc_dz         = Marshal.load(Marshal.dump(@condition.spbc_dz))

      # deep copy
      @condition = Marshal.load(Marshal.dump(c))
    end
  end

  class XYZ
    attr_accessor :temperature, :atoms, :size
    attr_reader :index_of_frame, :eof

    def initialize
      @index_of_frame = -1
      @eof            = false
    end

    # 1 frame
    def read fp
      @size        = fp.gets.to_i
      @eof         = true if @size == 0
      @temperature = fp.gets.to_f
      @atoms       = []
      @size.times do
        words           = fp.gets.split
        atom            = XYZAtom.new
        atom.type       = words[0]
        atom.coordinate = Coordinate.new(words[1].to_f, words[2].to_f, words[3].to_f)
        atom.energy     = words[4].to_f
        @atoms.push atom
      end
      @index_of_frame += 1
    end

    def write fp
      fp.puts @atoms.size, "%20.15e"%[@temperature]
      @atoms.each_with_index do |atom, index|
        fp.puts "#{atom.type} % 20.15e % 20.15e % 20.15e % 20.15e #{index}"%[atom.coordinate.x, atom.coordinate.y, atom.coordinate.z, atom.energy]
      end
    end
  end

  class Element
    attr_accessor :atom_type, :mass
  end

  class Condition
    attr_accessor :loop_number,:quench_steps, :annealing_steps,
      :dynamics_steps, :output_interval, :cg_steps, :dt,
      :temperature, :dtemperature, :elements, :potential_type,
      :displacement, :displacement_u, :displacement_ux, :displacement_uz, :height,
      :shear_stress, :stress, :stress_ex, :stress_ez,
      :spbc_dz, :periodic

    def read fp
      # 全体ループの数
      @loop_number     = fp.gets.to_i
      # 数値急冷法のステップ数
      @quench_steps    = fp.gets.to_i
      # 焼き鈍しステップ数（NVTアンサンブル）
      @annealing_steps = fp.gets.to_i
      # NVEアンサンブルステップ数
      nve              = fp.gets.split
      @dynamics_steps  = nve[0].to_i
      @output_interval = nve[1].to_i
      # 共役勾配法ステップ数
      @cg_steps        = fp.gets.to_i
      # 時間刻みの大きさ（単位：秒）
      @dt              = fp.gets.to_f
      # 設定温度と温度の変化量（焼き鈍しで使う）（単位：K）
      temperatures     = fp.gets.split
      @temperature     = temperatures[0].to_f
      @dtemperature    = temperatures[1].to_f
      # 計算に含まれる原子の種類数
      elements_num     = fp.gets.to_i
      @elements        = []
      elements_num.times do
        el             = Element.new
        # 原子のラベル（元素記号）と質量（単位：kg）
        words          = fp.gets.split
        el.atom_type   = words[0].chomp
        el.mass        = words[1].to_f
        @elements.push el
      end
      # ポテンシャルタイプの読み込み
      @potential_type  = fp.gets.chomp
      @potential_files = []
      case @potential_type
      when 'pair'
        summation(elements_num).times{ @potential_files.push fp.gets.chomp }
      when 'eam'
        (summation(elements_num) + elements_num * 2).times{ @potential_files.push fp.gets.chomp }
      when 'tersoff'
        summation(elements_num).times{ @potential_files.push fp.gets.chomp }
      when 'adp'
        (summation(elements_num) * 3 + elements_num * 2).times{ @potential_files.push fp.gets.chomp }
      when 'sw'
        @potential_files.push fp.gets.chomp
      when 'brenner'
        summation(elements_num).times{ @potential_files.push fp.gets.chomp }
      when '2bm'
        elements_num.times{ @potential_files.push fp.gets.chomp }
      else
        puts '未定義のポテンシャルタイプです。'
        exit false
      end
      # 変位境界条件（条件を与えるか、その大きさ、与える領域の厚さ）
      disp_bc          = fp.gets.split
      @displacement    = disp_bc[0].to_i.to_bool
      @displacement_u  = disp_bc[1].to_f
      @displacement_ux = disp_bc[2].to_f
      @displacement_uz = disp_bc[3].to_f
      @height          = disp_bc[4].to_f
      # 応力境界条件（条件を与えるか、その大きさ）
      shear_bc         = fp.gets.split
      @shear_stress    = shear_bc[0].to_i.to_bool
      @stress          = shear_bc[1].to_f
      @stress_ex       = shear_bc[2].to_f
      @stress_ez       = shear_bc[3].to_f
      # SPBC_dx
      @spbc_dz         = fp.gets.to_f
      # 周期境界条件を与えるか
      periodic_bc      = fp.gets.split
      @periodic        = Periodic.new(
        periodic_bc[0].to_i.to_bool,
        periodic_bc[1].to_i.to_bool,
        periodic_bc[2].to_i.to_bool
      )
    end

    def read_from_mdl fp
      # 全体ループの数
      @loop_number     = fp.gets.to_i
      # 数値急冷法のステップ数
      @quench_steps    = fp.gets.to_i
      # 焼き鈍しステップ数（NVTアンサンブル）
      @annealing_steps = fp.gets.to_i
      # NVEアンサンブルステップ数
      nve = fp.gets.split
      @dynamics_steps  = nve[0].to_i
      @output_interval = nve[1].to_i
      # 共役勾配法ステップ数
      @cg_steps        = fp.gets.to_i
      # 時間刻みの大きさ（単位：秒）
      @dt              = fp.gets.to_f
      # 設定温度と温度の変化量（焼き鈍しで使う）（単位：K）
      temperatures     = fp.gets.split
      @temperature     = temperatures[0].to_f
      @dtemperature    = temperatures[1].to_f
      # 計算に含まれる原子の種類数
      elements_num     = fp.gets.to_i
      @elements        = []
      elements_num.times do
        el             = Element.new
        # 原子のラベル（元素記号）と質量（単位：kg）
        words          = fp.gets.split
        el.atom_type   = words[0].chomp
        el.mass        = words[1].to_f
        @elements.push el
      end
      # ポテンシャルタイプの読み込み
      @potential_type  = fp.gets.chomp
      @potential_files = []
      case @potential_type
      when 'pair'
        summation(elements_num).times{ @potential_files.push fp.gets.chomp }
      when 'eam'
        (summation(elements_num) + elements_num * 2).times{ @potential_files.push fp.gets.chomp }
      when 'tersoff'
        summation(elements_num).times{ @potential_files.push fp.gets.chomp }
      when 'adp'
        (summation(elements_num) * 3 + elements_num * 2).times{ @potential_files.push fp.gets.chomp }
      when 'sw'
        @potential_files.push fp.gets.chomp
      when 'brenner'
        summation(elements_num).times{ @potential_files.push fp.gets.chomp }
      when '2bm'
        elements_num.times{ @potential_files.push fp.gets.chomp }
      else
        puts '未定義のポテンシャルタイプです。'
        exit false
      end
      # 原子の数
      size             = fp.gets.to_i
      # 計算領域の大きさ
      volumes          = fp.gets.split
      volume           = Volume.new(
        volumes[0].to_f,
        volumes[1].to_f,
        volumes[2].to_f
      )
      # 変位境界条件（条件を与えるか、その大きさ、与える領域の厚さ）
      disp_bc          = fp.gets.split
      @displacement    = disp_bc[0].to_i.to_bool
      @displacement_u  = disp_bc[1].to_f
      @displacement_ux = disp_bc[2].to_f
      @displacement_uz = disp_bc[3].to_f
      @height          = disp_bc[4].to_f
      # 応力境界条件（条件を与えるか、その大きさ）
      shear_bc         = fp.gets.split
      @shear_stress    = shear_bc[0].to_i.to_bool
      @stress          = shear_bc[1].to_f
      @stress_ex       = shear_bc[2].to_f
      @stress_ez       = shear_bc[3].to_f
      # SPBC_dx
      @spbc_dz         = fp.gets.to_f
      # 周期境界条件を与えるか
      periodic_bc      = fp.gets.split
      @periodic        = Periodic.new(
        periodic_bc[0].to_i.to_bool,
        periodic_bc[1].to_i.to_bool,
        periodic_bc[2].to_i.to_bool
      )
      return volume, size
    end

    def write fp
      fp.puts @loop_number,
        @quench_steps,
        @annealing_steps,
        "#{@dynamics_steps} #{@output_interval}",
        @cg_steps,
        "%20.15e"%[@dt],
        "%20.15e %20.15e"%[@temperature, @dtemperature],
        @elements.size
      @elements.each { |el| fp.puts "#{el.atom_type} %20.15e"%[el.mass] }
      fp.puts @potential_type
      @potential_files.each { |pf| fp.puts pf }
      fp.puts "#{@displacement.to_i} %20.15e %20.15e %20.15e %20.15e"%[@displacement_u, @displacement_ux, @displacement_uz, @height],
        "#{@shear_stress.to_i} %20.15e %20.15e %20.15e"%[@stress, @stress_ex, @stress_ez],
        "%20.15e"%[@spbc_dz],
        "#{@periodic.x.to_i} #{@periodic.y.to_i} #{@periodic.z.to_i}"
    end

    def write_to_mdl fp, volume, size
      fp.puts @loop_number,
        @quench_steps,
        @annealing_steps,
        "#{@dynamics_steps} #{@output_interval}",
        @cg_steps,
        "%20.15e"%[@dt],
        "%20.15e %20.15e"%[@temperature, @dtemperature],
        @elements.size
      @elements.each { |el| fp.puts "#{el.atom_type} %20.15e"%[el.mass] }
      fp.puts @potential_type
      @potential_files.each { |pf| fp.puts pf }
      fp.puts size,
        "%20.15e %20.15e %20.15e"%[volume.x, volume.y, volume.z],
        "#{@displacement.to_i} %20.15e %20.15e %20.15e %20.15e"%[@displacement_u, @displacement_ux, @displacement_uz, @height],
        "#{@shear_stress.to_i} %20.15e %20.15e %20.15e"%[@stress, @stress_ex, @stress_ez],
        "%20.15e"%[@spbc_dz],
        "#{@periodic.x.to_i} #{@periodic.y.to_i} #{@periodic.z.to_i}"
    end
  end
end
